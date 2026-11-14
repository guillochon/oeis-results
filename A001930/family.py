"""Every OEIS entry that follows from the cycle index series of posets (degree <= 17) or from the labeled
counts A001035 (n <= 19).  Checks each against the OEIS terms in KNOWN, prints the new terms, and writes
family.json, family.md and the b-files.

    python family.py
"""
import json, os, subprocess, sys
from fractions import Fraction
from math import factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = [sys.argv[0], "17"]
from topologies import load_fix, partitions, fix, z, count, inv_euler, A000112, A001035  # noqa: E402
import prime as sf  # noqa: E402

N = 17
P19 = A001035 + [646099441937791106493755218560442089979]  # A001035(19), Ayala 2026 (in the b-file)

# OEIS terms as of the Oct 23 2026 snapshot: (offset, terms)
KNOWN = {
    "A001930": (0, [1, 1, 3, 9, 33, 139, 718, 4535, 35979, 363083, 4717687, 79501654, 1744252509, 49872339897, 1856792610995, 89847422244493, 5637294117525695]),
    "A001928": (0, [1, 1, 2, 6, 21, 94, 512, 3485, 29515, 314474, 4255727, 73831813, 1653083021, 47941962135, 1803010446411, 87882300251730, 5543501326580737]),
    "A326898": (0, [1, 2, 5, 14, 47, 186, 904, 5439, 41418, 404501, 5122188, 84623842, 1828876351, 51701216248, 1908493827243, 91755916071736, 5729050033597431]),
    "A091073": (0, [1, 2, 8, 39, 242, 1895, 19051, 246895, 4145108, 90325655, 2555630036, 93810648902, 4461086120602, 274339212258846, 21775814889230580, 2226876304576948549]),
    "A296105": (0, [1, 2, 5, 25, 157, 1325, 14358, 199763, 3549001, 80673244, 2352747542, 88240542454, 4261209044877, 264988507673267, 21207485269909946, 2182146922863398203]),
    "A399304": (0, [1, 1, 6, 31, 203, 1653, 17156, 227844, 3898213, 86180547, 2465304381, 91255018866, 4367275471700, 269878126138244, 21501475676971734, 2205100489687717969]),
    "A079265": (0, [1, 2, 7, 32, 192, 1490, 15067, 198296, 3398105, 75734592, 2191591226, 82178300654, 3984499220967, 249298391641352, 20089200308020179, 2081351202770089728]),
    "A046907": (0, [1, 1, 1, 2, 7, 31, 184, 1351, 12524, 146468, 2177570, 41374407, 1008220289, 31559446774, 1269310589336, 65562045668340, 4345161435996517]),
    "A046911": (0, [1, 1, 2, 4, 14, 62, 373, 2722]),
    "A046909": (0, [1, 1, 1, 1, 2, 17, 167, 1672]),
    "A046904": (0, [1, 1, 0, 0, 1, 4, 28, 234]),
    "A354693": (1, [1, 0, 0, 1, 4, 28, 234, 2585, 36326, 646405, 14528011, 412212506]),
    "A066304": (0, [1, 1, 1, 3, 8, 30, 150]),
    "A066305": (1, [1, 0, 2, 5, 22, 117]),
    "A381121": (0, [1, 0, 1, 3, 11, 47, 255, 1727, 14954, 166232, 2384053, 44182143, 1058142319, 32718935706, 1304369332319, 66936884741385, 4414855587293931]),
    "A028849": (2, [0, 0, 6, 0, 27, 0, 91, 45, 417, 0]),
    "A028850": (2, [0, 0, 3, 0, 10, 0, 30, 15, 126, 0]),
    "A028851": (2, [0, 0, 3, 0, 12, 0, 40, 21, 188, 0]),
    "A028852": (2, [0, 0, 1, 0, 3, 0, 10, 6, 44, 0]),
    "A028855": (2, [3, 9, 27, 139, 691, 4535]),
    "A028856": (2, [2, 5, 13, 63, 308, 2045, 16969, 183216, 2567158, 46749427]),
    "A028857": (2, [2, 6, 18, 94, 500, 3485]),
    "A028858": (2, [1, 3, 9, 44, 235, 1650, 14502, 163335, 2360675, 43944974]),
    # labeled
    "A046905": (0, [1, 1, 0, 0, 24, 360, 17400, 1066800, 97463520]),  # a(8) added to the OEIS Nov 7 2026
    "A066302": (0, [1, 1, 1, 7, 75, 1531, 50673, 2613703, 202180723, 22853355895, 3701983130913, 846741042881779, 270316009546766571, 119343586350194910211, 72325998629777739416209, 59798727157327673157936319]),
    "A066303": (1, [1, 0, 6, 50, 1220, 42122, 2278248, 182111666, 21094774212, 3479970392642, 807001730170592, 260359119927269882, 115887302333840512956, 70677462143507496419690, 58725460618979611312258632]),
    "A085628": (0, [1, 2, 12, 152, 3504, 135392, 8321472, 784621952, 110521185024, 22789653765632, 6769730814753792, 2859584874712881152, 1699286839524775931904, 1407801166901961190203392, 1613567168628788544015286272, 2541721059997800475952740401152, 5470980000021882982488097199161344, 16010777886758708558291439618510159872, 63422960133000705239390046164899600269312]),
    "A046906": (0, [1, 1, 0, 0, 24, 1080, 52440, 3281880, 277953144, 32418855000, 5239070305080, 1173944480658840, 363936227764858584, 155521768202208047640, 91218870039317505477720, 73113879800794757415243480, 79743817918540500914682249144, 117883366412734188786535902826200, 235329353612778837110901775412557560]),
    "A046908": (0, [1, 1, 1, 7, 97, 2251, 80821, 4305127, 332273257, 36630174931, 5711638291981, 1249898984911567, 381230073532620577, 161042140788424003291, 93667063572594041040421, 74610767840852891620692727, 80997478506602342803118178457, 119313601058907927882431190269731, 237541348427311374857037021264415741]),
    "A001927": (0, [1, 1, 2, 12, 146, 3060, 101642, 5106612, 377403266, 40299722580, 6138497261882, 1320327172853172, 397571105288091506, 166330355795371103700, 96036130723851671469482, 76070282980382554147600692, 82226869197428315925408327266, 120722306604121583767045993825620, 239727397782668638856762574296226842]),
}
NAMES = {
    "A001930": "topologies (quasi-orders), unlabeled",
    "A001928": "connected topologies, unlabeled",
    "A326898": "unlabeled topologies with up to n points",
    "A091073": "transitive relations, unlabeled",
    "A296105": "connected transitive relations, unlabeled",
    "A399304": "transitive relations with no isolated point, unlabeled",
    "A079265": "antisymmetric transitive relations, unlabeled",
    "A046907": "irreducible posets (not an ordinal sum), unlabeled",
    "A046911": "irreducible quasi-orders, unlabeled",
    "A046909": "connected irreducible quasi-orders, unlabeled",
    "A046904": "prime posets (no nontrivial module), unlabeled",
    "A354693": "prime posets, unlabeled",
    "A066304": "reduced posets (no two-point chain module), unlabeled",
    "A066305": "connected reduced posets, unlabeled",
    "A381121": "posets covering n points (first differences of A000112)",
    "A028849": "composite topologies", "A028850": "composite T_0 topologies",
    "A028851": "composite connected topologies", "A028852": "composite connected T_0 topologies",
    "A028855": "prime topologies", "A028856": "prime T_0 topologies",
    "A028857": "prime connected topologies", "A028858": "prime connected T_0 topologies",
    "A046905": "prime posets, labeled",
    "A066302": "reduced posets, labeled",
    "A066303": "connected reduced posets, labeled",
    "A085628": "antisymmetric transitive relations, labeled",
    "A046906": "connected irreducible posets, labeled",
    "A046908": "irreducible posets, labeled",
    "A001927": "connected posets, labeled",
}


# ---------- ordinary power series helpers (integer / Fraction lists) ----------
def smul(a, b, n):
    c = [0] * (n + 1)
    for i, x in enumerate(a[: n + 1]):
        if x:
            for j in range(min(len(b), n + 1 - i)):
                c[i + j] += x * b[j]
    return c


def sinv(a, n):  # 1/a, a[0] = 1
    b = [Fraction(1)] + [Fraction(0)] * n
    for m in range(1, n + 1):
        b[m] = -sum(a[k] * b[m - k] for k in range(1, min(m, len(a) - 1) + 1))
    return b


def slog(a, n):  # log a, a[0] = 1
    d = [Fraction(k * a[k]) for k in range(1, n + 1)]  # a'
    q = smul(d, sinv(a, n), n - 1)
    return [Fraction(0)] + [q[k - 1] / k for k in range(1, n + 1)]


def scompose(f, g, n):  # f(g(x)), g[0] = 0
    out, pw = [Fraction(0)] * (n + 1), [Fraction(1)] + [Fraction(0)] * n
    for k in range(n + 1):
        if k < len(f) and f[k]:
            for i in range(n + 1):
                out[i] += f[k] * pw[i]
        pw = smul(pw, g, n)
    return out


def srevert(g, n):  # h with g(h(x)) = x, g = x + ...
    h = [Fraction(0), Fraction(1)] + [Fraction(0)] * (n - 1)
    for m in range(2, n + 1):
        c = scompose(g, h[: m] + [Fraction(0)] * (n + 1 - m), m)
        h[m] = -c[m]
    return h


def ints(v):
    assert all(Fraction(x).denominator == 1 for x in v), v
    return [int(x) for x in v]


def dirichlet_primes(total, n):
    """prime(m), m = 2..n, when total(m) = #multisets of primes whose sizes multiply to m."""
    pr = {}
    ways = {1: 1}  # ways[m] = multisets using the primes found so far
    for m in range(2, n + 1):
        pr[m] = total[m] - ways.get(m, 0)
        assert pr[m] >= 0
        # multiply the Dirichlet series by (1 - m^-s)^(-pr[m])
        new = dict(ways)
        k, c = 1, 1
        while m ** k <= n:
            c = c * (pr[m] + k - 1) // k  # multisets of size k
            for d, w in ways.items():
                if d * m ** k <= n:
                    new[d * m ** k] = new.get(d * m ** k, 0) + w * c
            k += 1
        ways = new
    return pr


def main():
    tab = load_fix()
    out = {}
    pos, _ = count("poset", N, tab)
    assert pos == A000112
    top, ttri = count("top", N, tab)
    tra, _ = count("trans", N, tab)
    ctop = [1] + inv_euler(top)[1:]
    ctra = [1] + inv_euler(tra)[1:]
    cpos = [1] + inv_euler(pos)[1:]
    out["A001930"] = (0, top)
    out["A001928"] = (0, ctop)
    out["A326898"] = (0, [sum(top[: i + 1]) for i in range(N + 1)])
    out["A091073"] = (0, tra)
    out["A296105"] = (0, ctra)
    out["A399304"] = (0, [1] + [tra[i] - tra[i - 1] for i in range(1, N + 1)])
    out["A381121"] = (0, [1] + [pos[i] - pos[i - 1] for i in range(1, N + 1)])
    # antisymmetric transitive relations = posets with each point reflexive or not: Z_P[2x]
    ZP = {lam: Fraction(fix(lam, tab), z(lam)) for k in range(1, N + 1) for lam in partitions(k)}
    a = [Fraction(0)] * (N + 1)
    a[0] = Fraction(1)
    for lam, v in ZP.items():
        a[sum(lam)] += v * 2 ** len(lam)
    out["A079265"] = (0, ints(a))
    # irreducible = not an ordinal sum: g.f. 2 - 1/B(x)
    irr = lambda b: ints([1] + [-x for x in sinv(b, N)[1:]])
    out["A046907"] = (0, irr(pos))
    iq = irr(top)
    out["A046911"] = (0, iq)
    # a disconnected quasi-order on >= 2 points is irreducible
    out["A046909"] = (0, [1, 1] + [iq[n] - (top[n] - ctop[n]) for n in range(2, N + 1)])
    # prime posets: Z_Pr = (Z_C + Z_S - Z_P - p_1) o Z_P^<-1>
    ZC, ZS = sf.connected(ZP), sf.geometric_inverse(ZP)
    assert sf.count_unlabeled(ZC)[1:] == cpos[1:] and sf.count_unlabeled(ZS)[1:] == out["A046907"][1][1:]
    T = sf.add(sf.add(sf.add(ZC, ZS), ZP, -1), {(1,): Fraction(1)}, -1)
    ZPr = sf.right_divide(T, ZP, lo=4)
    pr = sf.count_unlabeled(ZPr)
    out["A046904"] = (0, [1, 1] + pr[2:])
    out["A354693"] = (1, [1] + pr[2:])
    # reduced posets: P = R o L_+ (chains substituted for points), Z_{L_+} = p_1/(1-p_1)
    ZR = sf.right_divide(ZP, {tuple([1] * m): Fraction(1) for m in range(1, N + 1)})
    red = sf.count_unlabeled(ZR)
    out["A066304"] = (0, [1] + red[1:])
    out["A066305"] = (1, inv_euler([1] + red[1:])[1:])
    # Kreminski: prime = not a product of two spaces with >= 2 points each
    for tot, pa, ca in ((top, "A028855", "A028849"), (pos, "A028856", "A028850"),
                        (ctop, "A028857", "A028851"), (cpos, "A028858", "A028852")):
        p = dirichlet_primes(tot, N)
        out[pa] = (2, [p[m] for m in range(2, N + 1)])
        out[ca] = (2, [tot[m] - p[m] for m in range(2, N + 1)])

    # ---------- labeled, n <= 19, from A001035 alone ----------
    M = 19
    B = [Fraction(P19[n], factorial(n)) for n in range(M + 1)]
    egf = lambda s: ints([s[n] * factorial(n) for n in range(M + 1)])
    Bp = [Fraction(0)] + B[1:]
    C = slog(B, M)
    S = [Fraction(0)] + [-x for x in sinv(B, M)[1:]]  # 1 - 1/B
    Tl = [C[n] + S[n] - Bp[n] - (1 if n == 1 else 0) for n in range(M + 1)]
    Prl = scompose(Tl, srevert(Bp, M), M)
    lp = egf(Prl)
    assert lp[:17] == sf.count_labeled(ZPr)[:17] + [0] * 0 or lp[: N + 1] == sf.count_labeled(ZPr)
    out["A046905"] = (0, [1, 1] + lp[2:])
    xo = [Fraction(0)] + [Fraction((-1) ** (k + 1)) for k in range(1, M + 1)]  # x/(1+x)
    R = scompose(B, xo, M)
    out["A066302"] = (0, egf(R))
    assert egf(R)[1 : N + 1] == sf.count_labeled(ZR)[1:]
    out["A066303"] = (1, egf(slog(R, M))[1:])
    out["A085628"] = (0, [2 ** n * P19[n] for n in range(M + 1)])
    out["A001927"] = (0, [1] + egf(C)[1:])
    out["A046908"] = (0, [1] + egf(S)[1:])
    # connected irreducible, labeled: irreducible minus disconnected (n >= 2)
    out["A046906"] = (0, [1, 1] + [egf(S)[n] - (P19[n] - egf(C)[n]) for n in range(2, M + 1)])

    # ---------- compare ----------
    rows, ok = [], True
    for a, (off, v) in out.items():
        ko, kv = KNOWN[a]
        assert ko == off, a
        m = min(len(kv), len(v))
        match = kv[:m] == v[:m]
        ok &= match
        new = v[len(kv):]
        rows.append((a, off, len(kv), match, new))
        print(f"{a} {NAMES[a]}: {'OK' if match else 'MISMATCH'} on {m} known terms; new: "
              + (", ".join(f"a({off + len(kv) + i}) = {x}" for i, x in enumerate(new)) or "none"))
        if not match:
            print("   ours:", v[:m], "\n   OEIS:", kv[:m])
    json.dump({a: {"offset": off, "terms": v} for a, (off, v) in out.items()},
              open(os.path.join(HERE, "family.json"), "w"), indent=0)
    json.dump({"N": N, "terms": [[list(k), str(v)] for k, v in sorted(ZP.items(), key=lambda kv: (sum(kv[0]), kv[0]))]},
              open(os.path.join(HERE, "poset_cycle_index.json"), "w"))
    with open(os.path.join(HERE, "family.md"), "w", encoding="utf-8") as fh:
        fh.write("| Entry | What | Known terms reproduced | New terms |\n|---|---|---|---|\n")
        for a, off, nk, match, new in rows:
            fh.write(f"| {a} | {NAMES[a]} | {nk} ({'all' if match else 'MISMATCH'}) | "
                     + ("<br>".join(f"a({off + nk + i}) = {x}" for i, x in enumerate(new)) or "none") + " |\n")
    for a, (off, v) in out.items():
        if len(v) > 15 and len(v) > len(KNOWN[a][1]):
            with open(os.path.join(HERE, f"b{a[1:]}.txt"), "w", newline="\n") as fh:
                for i, x in enumerate(v):
                    fh.write(f"{off + i} {x}\n")
    print("ALL KNOWN TERMS REPRODUCED" if ok else "SOME MISMATCH")


if __name__ == "__main__":
    main()
