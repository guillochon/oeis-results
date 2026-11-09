"""Beam search for characters chi: primes -> Z/k with a late first zero pair.

chi(q) = 1 for every prime q except a small exception set E (q -> v). chi(n) = sum_q chi(q) nu_q(n) mod k.
The score of E is the first r with chi(r) = chi(r+1) = 0 (so D(k, r) holds and a(k) >= r).
Moves: change chi(q) for a prime q dividing r or r+1 at the current first pair.

usage: python ls/beam.py k X width rounds start_exc... (e.g. 3:0)  -> writes ls/best_k{k}.txt
"""
import sys, time
import numpy as np

k, X, W, R = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
start = {int(a.split(':')[0]): int(a.split(':')[1]) % k for a in sys.argv[5:]}
t0 = time.time()
M = X + 2
rem = np.arange(M, dtype=np.int32)
chi = np.zeros(M, dtype=np.int8)
lim = int(M ** 0.5) + 1
sv = np.ones(lim + 1, bool); sv[:2] = False
for i in range(2, int(lim ** 0.5) + 1):
    if sv[i]: sv[i*i::i] = False
small = [int(p) for p in np.nonzero(sv)[0]]
for p in small:
    pe = p
    while pe < M:
        rem[pe::pe] //= p; chi[pe::pe] += 1; pe *= p
chi[2:][rem[2:] > 1] += 1
del rem
chi %= k
print(f"base built {time.time()-t0:.1f}s", flush=True)

def add(p, d):
    if d == 0: return
    pe = p
    while pe < M:
        s = chi[pe::pe]; s += d; s %= k   # int8 view, values stay < 2k
        pe *= p

def first_pair():
    z = chi == 0
    z[0] = False
    w = z[1:-1] & z[2:]
    i = int(np.argmax(w))
    return i + 1 if w[i] else M

def factor(n):
    f = []; d = 2
    while d * d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: f.append(n)
    return f

def apply(E, sign):
    for p, v in E.items():
        add(p, (sign * (v - 1)) % k)

def score(E):
    apply(E, 1); r = first_pair(); apply(E, -1); return r

beam = [(score(start), start)]
seen = {frozenset(start.items())}
best = beam[0]
print(f"start {best}", flush=True)
for rnd in range(R):
    cand = []
    for r, E in beam:
        apply(E, 1)
        for p in factor(r) + factor(r + 1):
            cur = E.get(p, 1)
            for v in range(k):
                if v == cur: continue
                E2 = dict(E)
                if v == 1: E2.pop(p)
                else: E2[p] = v
                key = frozenset(E2.items())
                if key in seen: continue
                seen.add(key)
                add(p, (v - cur) % k)
                cand.append((first_pair(), E2))
                add(p, (cur - v) % k)
        apply(E, -1)
    if not cand: break
    cand.sort(key=lambda t: (-t[0], len(t[1])))
    beam = cand[:W]
    if beam[0][0] > best[0]:
        best = beam[0]
        with open(f"ls/best_k{k}.txt", "w") as f:
            f.write(f"{best[0]} {sorted(best[1].items())}\n")
    print(f"round {rnd} {time.time()-t0:.0f}s top={[c[0] for c in beam[:5]]} best={best[0]} E={sorted(best[1].items())}", flush=True)
    if best[0] >= M: break
