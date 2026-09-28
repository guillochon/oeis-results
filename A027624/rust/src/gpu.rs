//! Native GPU driver for a(8): runs the CUDA kernel of `../cuda/a8.cu` (embedded as a cubin) through
//! the CUDA driver API in `nvcuda.dll`, loaded at run time. No CUDA toolkit is needed on Windows,
//! only the NVIDIA driver; the cubin is rebuilt in WSL with `cuda/ptx.sh` when the kernel changes.
//!
//! This is a port of the host half of `a8.cu`: the same tables, batches, output format (the result
//! files are interchangeable with the WSL program), stall watchdog, safe resume and check mode.
//!
//!   hypercube-indep a8 [--from I] [--to J] [--stride S] [--batch B] [--sym] [--out FILE]
//!                      [--check REF.txt] [--stall SECS] [--reps FILE]
//!   hypercube-indep a8run [--nosym] [--batch B] [--stall SECS]
//!       the full run: 4000-orbit chunks, each in a child `a8` process that is restarted if it
//!       stalls or times out; log in ../logs/a8_gpu.log

use crate::ygroup::{load_reps6, U256};
use std::collections::HashSet;
use std::ffi::c_void;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

static CUBIN: &[u8] = include_bytes!("../../cuda/a8.cubin");

const NI4: usize = 743;
const NI3: usize = 35;
const THREADS: u32 = 256; // must match THREADS in a8.cu

// ---------------------------------------------------------------------------------------------
// CUDA driver API, loaded from nvcuda.dll
// ---------------------------------------------------------------------------------------------

type CuResult = i32;
type DevPtr = u64;
type Handle = *mut c_void;
const CUDA_ERROR_NOT_READY: CuResult = 600;

#[allow(non_snake_case)]
struct Cuda {
    cuInit: unsafe extern "system" fn(u32) -> CuResult,
    cuDeviceGet: unsafe extern "system" fn(*mut i32, i32) -> CuResult,
    cuDevicePrimaryCtxRetain: unsafe extern "system" fn(*mut Handle, i32) -> CuResult,
    cuCtxSetCurrent: unsafe extern "system" fn(Handle) -> CuResult,
    cuModuleLoadData: unsafe extern "system" fn(*mut Handle, *const c_void) -> CuResult,
    cuModuleGetFunction: unsafe extern "system" fn(*mut Handle, Handle, *const u8) -> CuResult,
    cuModuleGetGlobal: unsafe extern "system" fn(*mut DevPtr, *mut usize, Handle, *const u8) -> CuResult,
    cuMemAlloc: unsafe extern "system" fn(*mut DevPtr, usize) -> CuResult,
    cuMemcpyHtoD: unsafe extern "system" fn(DevPtr, *const c_void, usize) -> CuResult,
    cuMemcpyDtoH: unsafe extern "system" fn(*mut c_void, DevPtr, usize) -> CuResult,
    cuLaunchKernel: unsafe extern "system" fn(
        Handle, u32, u32, u32, u32, u32, u32, u32, Handle, *mut *mut c_void, *mut *mut c_void,
    ) -> CuResult,
    cuEventCreate: unsafe extern "system" fn(*mut Handle, u32) -> CuResult,
    cuEventRecord: unsafe extern "system" fn(Handle, Handle) -> CuResult,
    cuEventQuery: unsafe extern "system" fn(Handle) -> CuResult,
    cuGetErrorString: unsafe extern "system" fn(CuResult, *mut *const u8) -> CuResult,
}

#[cfg(windows)]
fn load_symbol(lib: isize, name: &str) -> *const c_void {
    #[link(name = "kernel32")]
    extern "system" {
        fn GetProcAddress(h: isize, name: *const u8) -> *const c_void;
    }
    let cname = format!("{name}\0");
    let p = unsafe { GetProcAddress(lib, cname.as_ptr()) };
    assert!(!p.is_null(), "nvcuda.dll has no {name}");
    p
}

impl Cuda {
    #[cfg(windows)]
    fn load() -> Cuda {
        #[link(name = "kernel32")]
        extern "system" {
            fn LoadLibraryA(name: *const u8) -> isize;
        }
        let lib = unsafe { LoadLibraryA(b"nvcuda.dll\0".as_ptr()) };
        assert!(lib != 0, "could not load nvcuda.dll (is the NVIDIA driver installed?)");
        macro_rules! sym {
            ($n:literal) => {
                unsafe { std::mem::transmute(load_symbol(lib, $n)) }
            };
        }
        Cuda {
            cuInit: sym!("cuInit"),
            cuDeviceGet: sym!("cuDeviceGet"),
            cuDevicePrimaryCtxRetain: sym!("cuDevicePrimaryCtxRetain"),
            cuCtxSetCurrent: sym!("cuCtxSetCurrent"),
            cuModuleLoadData: sym!("cuModuleLoadData"),
            cuModuleGetFunction: sym!("cuModuleGetFunction"),
            cuModuleGetGlobal: sym!("cuModuleGetGlobal_v2"),
            cuMemAlloc: sym!("cuMemAlloc_v2"),
            cuMemcpyHtoD: sym!("cuMemcpyHtoD_v2"),
            cuMemcpyDtoH: sym!("cuMemcpyDtoH_v2"),
            cuLaunchKernel: sym!("cuLaunchKernel"),
            cuEventCreate: sym!("cuEventCreate"),
            cuEventRecord: sym!("cuEventRecord"),
            cuEventQuery: sym!("cuEventQuery"),
            cuGetErrorString: sym!("cuGetErrorString"),
        }
    }

    #[cfg(not(windows))]
    fn load() -> Cuda {
        panic!("the native GPU driver is Windows-only; on Linux/WSL use cuda/a8")
    }

    fn check(&self, r: CuResult, what: &str) {
        if r != 0 {
            let mut s: *const u8 = std::ptr::null();
            unsafe { (self.cuGetErrorString)(r, &mut s) };
            let msg = if s.is_null() {
                "?".to_string()
            } else {
                unsafe { std::ffi::CStr::from_ptr(s as *const std::ffi::c_char) }.to_string_lossy().into_owned()
            };
            eprintln!("CUDA error {r} ({msg}) in {what}");
            std::process::exit(2);
        }
    }

    fn upload<T>(&self, data: &[T]) -> DevPtr {
        let bytes = std::mem::size_of_val(data);
        let mut p: DevPtr = 0;
        self.check(unsafe { (self.cuMemAlloc)(&mut p, bytes.max(1)) }, "cuMemAlloc");
        self.check(unsafe { (self.cuMemcpyHtoD)(p, data.as_ptr() as *const c_void, bytes) }, "cuMemcpyHtoD");
        p
    }

    fn set_const<T>(&self, module: Handle, name: &str, data: &[T]) {
        let cname = format!("{name}\0");
        let (mut p, mut size): (DevPtr, usize) = (0, 0);
        self.check(unsafe { (self.cuModuleGetGlobal)(&mut p, &mut size, module, cname.as_ptr()) }, name);
        let bytes = std::mem::size_of_val(data);
        assert_eq!(size, bytes, "size of constant {name}");
        self.check(unsafe { (self.cuMemcpyHtoD)(p, data.as_ptr() as *const c_void, bytes) }, name);
    }
}

// ---------------------------------------------------------------------------------------------
// Host tables (as in a8.cu)
// ---------------------------------------------------------------------------------------------

struct Tables {
    fm: Vec<u16>,
    i4: Vec<u16>,
    hid: Vec<u8>,
    lid: Vec<u8>,
    lstoff: Vec<u32>,
    lst: Vec<u16>,
    p3lo: Vec<u64>,
    p3hi: Vec<u64>,
    n4: Vec<u16>,
    hmask: Vec<u64>,
    i3: Vec<u8>,
    i3idx: Vec<u8>,
    border: Vec<u8>,
    aorder: Vec<u8>,
}

fn tables() -> Tables {
    let n4: Vec<u16> = (0..16usize).map(|p| (0..4).fold(0u16, |a, j| a | 1 << (p ^ 1 << j))).collect();
    let nb4 = |y: u32| (0..16).filter(|p| y >> p & 1 == 1).fold(0u32, |a, p| a | n4[p] as u32);
    let i4: Vec<u16> = (0..65536u32).filter(|&s| s & nb4(s) == 0).map(|s| s as u16).collect();
    assert_eq!(i4.len(), NI4);
    let mut f = vec![0u32; 65536];
    for &s in &i4 {
        f[s as usize] = 1;
    }
    for b in 0..16 {
        for x in 0..65536usize {
            if x >> b & 1 == 1 {
                f[x] += f[x ^ 1 << b];
            }
        }
    }
    let fm: Vec<u16> = f.iter().map(|&x| x as u16).collect();
    let i3: Vec<u8> = (0..256u32)
        .filter(|&s| {
            let n = (0..8).filter(|p| s >> p & 1 == 1).fold(0u32, |a, p| a | (0..3).fold(0, |b, j| b | 1 << (p ^ 1 << j)));
            s & n == 0
        })
        .map(|s| s as u8)
        .collect();
    assert_eq!(i3.len(), NI3);
    let pos = |b: u8| i3.iter().position(|&x| x == b).expect("I_4 set not a pair of I_3 sets") as u8;
    let hid: Vec<u8> = i4.iter().map(|&t| pos((t >> 8) as u8)).collect();
    let lid: Vec<u8> = i4.iter().map(|&t| pos((t & 255) as u8)).collect();
    let hmask: Vec<u64> = (0..256u32)
        .map(|h| (0..NI3).filter(|&i| i3[i] as u32 & !h == 0).fold(0u64, |a, i| a | 1 << i))
        .collect();
    let mut lstoff = vec![0u32; 65537];
    let mut lst = vec![];
    for x in 0..65536u32 {
        lst.extend(i4.iter().filter(|&&s| s as u32 & !x == 0));
        lstoff[x as usize + 1] = lst.len() as u32;
    }
    let p3: Vec<u128> = (0..=64).map(|i| 3u128.pow(i)).collect();
    let mut i3idx = vec![255u8; 256];
    for (i, &b) in i3.iter().enumerate() {
        i3idx[b as usize] = i as u8;
    }
    let mut border: Vec<u8> = (0..=255u8).collect();
    border.sort_by_key(|b| b.count_ones()); // stable, as std::stable_sort in a8.cu
    let mut aorder: Vec<u8> = (0..=255u8).collect();
    aorder.sort_by_key(|b| std::cmp::Reverse(b.count_ones()));
    Tables {
        fm,
        i4,
        hid,
        lid,
        lstoff,
        lst,
        p3lo: p3.iter().map(|&v| v as u64).collect(),
        p3hi: p3.iter().map(|&v| (v >> 64) as u64).collect(),
        n4,
        hmask,
        i3,
        i3idx,
        border,
        aorder,
    }
}

// ---------------------------------------------------------------------------------------------
// One run over a list of rep indices
// ---------------------------------------------------------------------------------------------

fn manifest() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

/// Drops a partial last line (cut off by a kill) and returns the indices with complete lines.
fn resume_set(path: &Path) -> HashSet<usize> {
    let mut done = HashSet::new();
    let Ok(all) = std::fs::read(path) else { return done };
    let keep = all.iter().rposition(|&c| c == b'\n').map(|p| p + 1).unwrap_or(0);
    if keep != all.len() {
        eprintln!("dropping a partial last line ({} bytes) from {}", all.len() - keep, path.display());
        std::fs::OpenOptions::new().write(true).open(path).unwrap().set_len(keep as u64).unwrap();
    }
    for line in String::from_utf8_lossy(&all[..keep]).lines() {
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() == 4 && f[3].bytes().all(|c| c.is_ascii_digit()) {
            if let Ok(i) = f[0].parse() {
                done.insert(i);
            }
        }
    }
    done
}

pub fn a8(args: &[String]) -> bool {
    let opt = |name: &str| args.iter().position(|a| a == name).map(|i| args[i + 1].clone());
    let flag = |name: &str| args.iter().any(|a| a == name);
    let reps_path = opt("--reps").map(PathBuf::from).unwrap_or_else(|| manifest().join("results/ylo_reps_d6.bin"));
    let reps = load_reps6(&reps_path);
    let sym = flag("--sym");
    let batch: usize = opt("--batch").map(|s| s.parse().unwrap()).unwrap_or(4);
    let stall: f64 = opt("--stall").map(|s| s.parse().unwrap()).unwrap_or(120.0);
    let check = opt("--check");
    let out_path = opt("--out").map(PathBuf::from);
    eprintln!("{} reps loaded; swap symmetry {}", reps.len(), if sym { "on (weighted rep sums)" } else { "off (plain rep sums)" });

    // which reps
    let mut todo: Vec<usize> = vec![];
    let mut expect: Vec<String> = vec![];
    if let Some(c) = &check {
        for line in std::fs::read_to_string(c).expect("reference file").lines() {
            let f: Vec<&str> = line.split_whitespace().collect();
            if f.len() == 4 {
                todo.push(f[0].parse().unwrap());
                expect.push(f[3].to_string());
            }
        }
    } else {
        let from: usize = opt("--from").map(|s| s.parse().unwrap()).unwrap_or(0);
        let to: usize = opt("--to").map(|s| s.parse().unwrap()).unwrap_or(reps.len()).min(reps.len());
        let stride: usize = opt("--stride").map(|s| s.parse().unwrap()).unwrap_or(1);
        let done = out_path.as_deref().map(resume_set).unwrap_or_default();
        todo = (from..to).step_by(stride).filter(|i| !done.contains(i)).collect();
        eprintln!("{} reps to do ({} already in the output)", todo.len(), done.len());
    }
    if todo.is_empty() {
        return true;
    }

    // CUDA setup
    let t = tables();
    let cu = Cuda::load();
    let (mut dev, mut ctx, mut module, mut func, mut ev): (i32, Handle, Handle, Handle, Handle) =
        (0, std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut(), std::ptr::null_mut());
    unsafe {
        cu.check((cu.cuInit)(0), "cuInit");
        cu.check((cu.cuDeviceGet)(&mut dev, 0), "cuDeviceGet");
        cu.check((cu.cuDevicePrimaryCtxRetain)(&mut ctx, dev), "cuDevicePrimaryCtxRetain");
        cu.check((cu.cuCtxSetCurrent)(ctx), "cuCtxSetCurrent");
        cu.check((cu.cuModuleLoadData)(&mut module, CUBIN.as_ptr() as *const c_void), "cuModuleLoadData");
        cu.check((cu.cuModuleGetFunction)(&mut func, module, b"slice_kernel\0".as_ptr()), "cuModuleGetFunction");
        cu.check((cu.cuEventCreate)(&mut ev, 2 /* CU_EVENT_DISABLE_TIMING */), "cuEventCreate");
    }
    cu.set_const(module, "c_p3lo", &t.p3lo);
    cu.set_const(module, "c_p3hi", &t.p3hi);
    cu.set_const(module, "c_n4", &t.n4);
    cu.set_const(module, "c_hmask", &t.hmask);
    cu.set_const(module, "c_i3", &t.i3);
    cu.set_const(module, "c_i3idx", &t.i3idx);
    cu.set_const(module, "c_border", &t.border);
    cu.set_const(module, "c_aorder", &t.aorder);
    let d_fm = cu.upload(&t.fm);
    let d_im = cu.upload(&t.i4);
    let d_hid = cu.upload(&t.hid);
    let d_lid = cu.upload(&t.lid);
    let d_lstoff = cu.upload(&t.lstoff);
    let d_lst = cu.upload(&t.lst);
    let d_brep = cu.upload(&vec![0u32; batch]);
    let d_out = cu.upload(&vec![0u64; batch * 65536 * 3]);
    let mut h_out = vec![0u64; batch * 65536 * 3];

    let mut fo = out_path.as_ref().map(|p| std::fs::OpenOptions::new().create(true).append(true).open(p).unwrap());
    let mut bad = 0;
    let t0 = Instant::now();
    for (k0, chunk) in todo.chunks(batch).enumerate() {
        let k = k0 * batch;
        let nb = chunk.len();
        let br: Vec<u32> = chunk.iter().map(|&i| reps[i].0).collect();
        cu.check(unsafe { (cu.cuMemcpyHtoD)(d_brep, br.as_ptr() as *const c_void, nb * 4) }, "cuMemcpyHtoD reps");
        // kernel arguments, in the order of slice_kernel's parameters
        let (mut a_reps, mut a_rep0, mut a_fm, mut a_im, mut a_hid, mut a_lid, mut a_lstoff, mut a_lst, mut a_out, mut a_sym) =
            (d_brep, 0i32, d_fm, d_im, d_hid, d_lid, d_lstoff, d_lst, d_out, sym as i32);
        let mut params: [*mut c_void; 10] = [
            &mut a_reps as *mut _ as *mut c_void,
            &mut a_rep0 as *mut _ as *mut c_void,
            &mut a_fm as *mut _ as *mut c_void,
            &mut a_im as *mut _ as *mut c_void,
            &mut a_hid as *mut _ as *mut c_void,
            &mut a_lid as *mut _ as *mut c_void,
            &mut a_lstoff as *mut _ as *mut c_void,
            &mut a_lst as *mut _ as *mut c_void,
            &mut a_out as *mut _ as *mut c_void,
            &mut a_sym as *mut _ as *mut c_void,
        ];
        unsafe {
            cu.check(
                (cu.cuLaunchKernel)(func, 65536, nb as u32, 1, THREADS, 1, 1, 0, std::ptr::null_mut(), params.as_mut_ptr(), std::ptr::null_mut()),
                "cuLaunchKernel",
            );
            cu.check((cu.cuEventRecord)(ev, std::ptr::null_mut()), "cuEventRecord");
        }
        // stall watchdog (see a8.cu)
        let tl = Instant::now();
        loop {
            let q = unsafe { (cu.cuEventQuery)(ev) };
            if q == 0 {
                break;
            }
            if q != CUDA_ERROR_NOT_READY {
                cu.check(q, "cuEventQuery");
            }
            if tl.elapsed().as_secs_f64() > stall {
                eprintln!("STALL: batch starting at rep {} not finished after {stall:.0} s; exiting for a restart", chunk[0]);
                std::process::exit(3);
            }
            std::thread::sleep(Duration::from_millis(2));
        }
        cu.check(unsafe { (cu.cuMemcpyDtoH)(h_out.as_mut_ptr() as *mut c_void, d_out, nb * 65536 * 3 * 8) }, "cuMemcpyDtoH");
        let mut lines = String::new();
        for (j, &i) in chunk.iter().enumerate() {
            let mut s = U256::default();
            for o in h_out[j * 65536 * 3..(j + 1) * 65536 * 3].chunks_exact(3) {
                s.add(U256 { hi: o[2] as u128, lo: (o[1] as u128) << 64 | o[0] as u128 });
            }
            let d = s.to_decimal();
            if check.is_some() {
                let ok = d == expect[k + j];
                bad += !ok as usize;
                eprintln!("rep #{i}: gpu {d}  cpu {}  {}", expect[k + j], if ok { "OK" } else { "MISMATCH" });
            } else {
                lines += &format!("{} {} {} {}\n", i, reps[i].0, reps[i].1, d);
            }
        }
        match (&mut fo, lines.is_empty()) {
            (Some(f), false) => {
                f.write_all(lines.as_bytes()).unwrap();
                f.flush().unwrap();
            }
            (None, false) => print!("{lines}"),
            _ => {}
        }
        let done = k + nb;
        if check.is_none() && (done % (batch * 50) == 0 || done == todo.len()) {
            let el = t0.elapsed().as_secs_f64();
            eprintln!(
                "[{el:8.0}s] {done}/{} reps  ({:.3} s/rep, ETA {:.1} h)",
                todo.len(),
                el / done as f64,
                el / done as f64 * (todo.len() - done) as f64 / 3600.0
            );
        }
    }
    let el = t0.elapsed().as_secs_f64();
    eprintln!(
        "done: {} reps in {el:.1}s ({:.3} s/rep){}",
        todo.len(),
        el / todo.len() as f64,
        if check.is_none() { "" } else if bad > 0 { "  CHECK FAILED" } else { "  ALL MATCH" }
    );
    bad == 0
}

// ---------------------------------------------------------------------------------------------
// Supervisor for the full run
// ---------------------------------------------------------------------------------------------

pub fn a8run(args: &[String]) -> bool {
    let opt = |name: &str| args.iter().position(|a| a == name).map(|i| args[i + 1].clone());
    let nosym = args.iter().any(|a| a == "--nosym");
    let batch = opt("--batch").unwrap_or_else(|| "4".into());
    let stall = opt("--stall").unwrap_or_else(|| "120".into());
    let base = manifest().join("..");
    let out = opt("--out").map(PathBuf::from).unwrap_or_else(|| base.join(if nosym { "results/a8_plain.txt" } else { "results/a8_sym.txt" }));
    let reference = base.join(if nosym { "logs/cpu6_reference.txt" } else { "logs/cpu6_reference_sym.txt" });
    let log_path = opt("--log").map(PathBuf::from).unwrap_or_else(|| base.join("logs/a8_gpu.log"));
    std::fs::create_dir_all(base.join("results")).unwrap();
    let exe = std::env::current_exe().unwrap();
    let sym: Vec<&str> = if nosym { vec![] } else { vec!["--sym"] };
    let mut log = std::fs::OpenOptions::new().create(true).append(true).open(&log_path).unwrap();
    let mut say = |msg: String| {
        println!("{msg}");
        writeln!(log, "{msg}").ok();
    };

    // guard: the embedded kernel must reproduce the CPU reference sums
    let ok = std::process::Command::new(&exe)
        .args(["a8", "--check", reference.to_str().unwrap(), "--batch", &batch])
        .args(&sym)
        .stderr(std::process::Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false);
    if !ok {
        say("reference check failed; not starting".into());
        return false;
    }
    say(format!("native run: output {}, log {}", out.display(), log_path.display()));
    // --to / --chunk / --attempts / --out / --log exist for testing the supervisor on a small range
    let n: usize = opt("--to").map(|s| s.parse().unwrap()).unwrap_or(1_228_158);
    let ch: usize = opt("--chunk").map(|s| s.parse().unwrap()).unwrap_or(4000);
    let attempts: usize = opt("--attempts").map(|s| s.parse().unwrap()).unwrap_or(50);
    let t0 = Instant::now();
    for s in (0..n).step_by(ch) {
        for attempt in 1..=attempts {
            let logf = std::fs::OpenOptions::new().create(true).append(true).open(&log_path).unwrap();
            let mut child = std::process::Command::new(&exe)
                .args(["a8", "--from", &s.to_string(), "--to", &(s + ch).min(n).to_string(), "--batch", &batch, "--stall", &stall])
                .args(["--out", out.to_str().unwrap()])
                .args(&sym)
                .stdout(logf.try_clone().unwrap())
                .stderr(logf)
                .spawn()
                .expect("spawn a8");
            let tc = Instant::now();
            let status = loop {
                if let Some(st) = child.try_wait().unwrap() {
                    break Some(st);
                }
                if tc.elapsed() > Duration::from_secs(3600) {
                    child.kill().ok();
                    child.wait().ok();
                    break None;
                }
                std::thread::sleep(Duration::from_millis(500));
            };
            match status {
                Some(st) if st.success() => break,
                Some(st) => say(format!("chunk {s}: attempt {attempt} exited with {st}; retrying")),
                None => say(format!("chunk {s}: attempt {attempt} timed out; retrying")),
            }
            if attempt == attempts {
                say(format!("chunk {s}: giving up after {attempts} attempts"));
                return false;
            }
        }
        say(format!("[{:.0} s] reps done through {} of {n}", t0.elapsed().as_secs_f64(), (s + ch).min(n)));
    }
    say("finished; run aggregate.py".into());
    true
}
