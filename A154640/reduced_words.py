"""A154639: words in 5 letters avoiding the factors of Graham et al. (GLMWY III) / Mallows:
(i) V V, (ii) V1..V2m (m >= 2) with V1 = V3, V_{2j} = V_{2j+3} (1 <= j <= m-2), V_{2m-2} = V_{2m}."""
import sys
from itertools import product

def bad_factor(f):
    L = len(f)
    if L % 2 or L < 4: return False
    m = L // 2
    V = (None,) + tuple(f)
    return (V[1] == V[3] and all(V[2*j] == V[2*j+3] for j in range(1, m - 1))
            and V[2*m-2] == V[2*m])

def ok(w):
    if any(w[i] == w[i+1] for i in range(len(w) - 1)): return False
    return not any(bad_factor(w[i:j]) for i in range(len(w)) for j in range(i + 4, len(w) + 1, 2))

def brute(n):
    return sum(ok(w) for w in product(range(5), repeat=n))

if __name__ == "__main__":
    n = int(sys.argv[1])
    print([brute(k) for k in range(n + 1)])


def automaton(n):
    """DP over states (last three letters, set of parities of active partial matches).
    A partial match of the family (ii) needs, at its next relative index r:
    r odd (>= 5): letter equal to the letter three back; r even (>= 4): any
    letter, and the factor is complete (forbidden) if it equals the letter two back."""
    from collections import Counter
    cur = Counter({((), frozenset()): 1})
    out = [1]
    for _ in range(n):
        nxt = Counter()
        for (last, act), c in cur.items():
            for a in range(5):
                if last and a == last[-1]:
                    continue
                dead = False
                act2 = set()
                for par in act:  # par = parity of next relative index r
                    if par == 1:
                        if a == last[-3]:
                            act2.add(0)
                    else:
                        if a == last[-2]:
                            dead = True
                        act2.add(1)
                if dead:
                    continue
                if len(last) >= 2 and a == last[-2]:
                    act2.add(0)  # x y x seen: next relative index 4 (even)
                nxt[((last + (a,))[-3:], frozenset(act2))] += c
        cur = nxt
        out.append(sum(cur.values()))
    return out
