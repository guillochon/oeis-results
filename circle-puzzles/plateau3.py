"""a(n) = pieces of GG_n on the plateau above sec(pi/n), sampled at r = sec + f*(sec - 1), f in F; the first radius is
face-traced and area-checked.
Runs on half the cores at below-normal priority (stereohedra/pool.py); compares with the period-4 cubic fits."""
import ctypes, math, os, sys, time
from twodisk import ExactPuzzle as Puzzle

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "stereohedra"))
from pool import kernel32, make_pool  # noqa: E402


def no_ecoqos():
    """Opt out of power throttling, or Windows keeps background jobs on the E-cores."""
    if os.name == "nt":
        k = kernel32()
        state = (ctypes.c_ulong * 3)(1, 1, 0)   # Version 1, ControlMask EXECUTION_SPEED, StateMask 0
        k.SetProcessInformation(k.GetCurrentProcess(), 4, ctypes.byref(state), ctypes.sizeof(state))

FIT = {0: (1, 19, -32, 40), 1: (1, 73, -145, 127), 2: (1, -1, 0, 4), 3: (1, 77, -93, 79)}


def work(n):
    s = 1 / math.cos(math.pi / n)
    t = time.time()
    vals = []
    ok = None
    for f in F:
        p = Puzzle(n, s + f * (s - 1), 10**7)
        p.close()
        if ok is None:   # face-tracing check on the first radius: one outer face, areas sum to the union
            P, nb, nu, darea, _ = p.check()
            ok = nb == P and nu == 1 and abs(darea) < 1e-9
            vals.append(int(P))
        else:
            vals.append(int(p.pieces_fast()))
        del p
    a, b, c, d = FIT[n % 4]
    fit = (a * n**3 + b * n**2 + c * n + d) // 8
    return n, vals, fit, ok, time.time() - t


F = (1.3, 1.0, 1.6)

if __name__ == "__main__":
    ns = [int(x) for x in sys.argv[1].split(",")] if "," in sys.argv[1] else range(int(sys.argv[1]), int(sys.argv[2]) + 1)
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else None
    with make_pool(nproc, initializer=no_ecoqos) as pool:
        for n, vals, fit, ok, dt in pool.imap(work, ns):
            flag = "fit" if len(set(vals)) == 1 and vals[0] == fit else "DIFF"
            print(n, vals, fit, flag, "check-ok" if ok else "CHECK-FAIL", f"{dt:.1f}s", flush=True)
