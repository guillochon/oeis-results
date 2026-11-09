"""Fill the table T(n,k) (symmetric n x n arrays, every 8-connected region of 1s has >= k cells) by antidiagonals
d = n + k, running minregion_tm2 inside WSL.  Resumable: results go to tools/table_runs.txt ("n k value").

    python3 tools/run_table.py DMAX WORKERS MEM_GB

k = 1 is 2^(octant cells) and k = 2 comes from b_pinhole_d4.txt, so only k >= 3 is run.  Each job is limited to
MEM_GB of address space; a job that runs out is logged in tools/table_fail.txt, and any (n, k) with a failed
(n, k') for k' <= k or a failed (n', k) for n' <= n with the same parity class of octant size is skipped.
"""
import os, resource, subprocess, sys, threading, time
from concurrent.futures import ThreadPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
OUT, FAIL = os.path.join(HERE, 'table_runs.txt'), os.path.join(HERE, 'table_fail.txt')
DMAX, WORKERS, MEM = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
done = {}
if os.path.exists(OUT):
    for l in open(OUT):
        n, k, v = l.split(); done[(int(n), int(k))] = v
failed = set()
if os.path.exists(FAIL):
    for l in open(FAIL):
        n, k, g = l.split(); failed.add((int(n), int(k), float(g)))
lock = threading.Lock()

def dominated(n, k):
    return any(g >= MEM and fn <= n and fk <= k for fn, fk, g in failed)

def limit():
    b = int(MEM * 2**30); resource.setrlimit(resource.RLIMIT_AS, (b, b)); os.nice(10)

def job(nk):
    n, k = nk
    with lock:
        if nk in done or dominated(n, k): return
    t = time.time()
    r = subprocess.run([os.path.join(HERE, 'minregion_tm2'), str(n), str(k)], capture_output=True, text=True, preexec_fn=limit)
    with lock:
        if r.returncode == 0:
            line = r.stdout.strip(); done[nk] = line.split()[2]
            open(OUT, 'a').write(line + '\n')
            print(f'{line[:60]}  {time.time() - t:.0f}s  {r.stderr.strip()}', flush=True)
        else:
            failed.add((n, k, MEM)); open(FAIL, 'a').write(f'{n} {k} {MEM}\n')
            print(f'FAILED n={n} k={k} rc={r.returncode} {r.stderr.strip()[:80]}', flush=True)

jobs = [(d - k, k) for d in range(4, DMAX + 1) for k in range(3, d)]
jobs = [j for j in jobs if 8 * j[1] <= 511]
with ThreadPoolExecutor(WORKERS) as ex:
    list(ex.map(job, jobs))
print('finished', flush=True)
