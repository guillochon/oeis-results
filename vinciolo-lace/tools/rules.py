"""Window rules: an a x b window (any of the 8 symmetries, frame = filled) is allowed iff it occurs in the
given panels."""
import numpy as np
def codes(A, a, b):
    pad = max(a, b) - 1
    P = np.pad(A.astype(np.int64), pad, constant_values=1); m1 = P.shape[0] - a + 1; m2 = P.shape[1] - b + 1
    code = np.zeros((m1, m2), np.int64)
    for x in range(a):
        for y in range(b): code = code * 2 + P[x:x + m1, y:y + m2]
    return code
def d4(A):
    imgs = [np.rot90(A, k) for k in range(4)]; return imgs + [x.T for x in imgs]
def table(panels, a, b):
    ok = np.zeros(1 << (a * b), bool)
    for A in panels:
        for B in d4(A): ok[codes(B, a, b).ravel()] = True
    return ok
def violations(A, a, b, ok):
    return int((~ok[codes(A, a, b)]).sum() + (~ok[codes(A.T, a, b)]).sum())
