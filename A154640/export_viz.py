"""Export spheres (center, radius, strategy-(a) and -(b) generation) for the three.js page.

Geometry: Mallows's start in augmented curvature-center coordinates (bbar, b, b*x):
three unit spheres centered at (2/sqrt3, 0, 0), (-1/sqrt3, +-1, 0) and the planes
z = +-1.  A sphere with integer vector v is v * C_AUG.

Strategy (a) generation of a sphere = length of the minimal coset representative:
v is e_j^T M_w in the Tits cone; reflecting at any negative coordinate (v -> v R_k)
shortens the representative by one, until v = e_j.
"""
import base64, json, math, struct, sys
from engine import strategy_b, reflect, IDENT, N, bend

s3 = math.sqrt(3)
C_AUG = [
    (1 / 3, 1, 2 / s3, 0, 0),
    (1 / 3, 1, -1 / s3, 1, 0),
    (1 / 3, 1, -1 / s3, -1, 0),
    (2, 0, 0, 0, 1),
    (2, 0, 0, 0, -1),
]


def gen_a(v):
    v = list(v); steps = 0
    while True:
        k = next((k for k in range(N) if v[k] < 0), None)
        if k is None:
            assert sorted(v) == [0, 0, 0, 0, 1], v
            return steps
        # v R_k: coordinate k -> -v_k, coordinate m != k -> v_m + v_k
        vk = v[k]
        v = [(-vk if m == k else v[m] + vk) for m in range(N)]
        steps += 1


def geom(v):
    a = [sum(v[i] * C_AUG[i][c] for i in range(N)) for c in range(5)]
    b = a[1]
    return (a[2] / b, a[3] / b, a[4] / b, 1 / b)


if __name__ == "__main__":
    G = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    _, _, gb = strategy_b(G)
    rows = []
    for v, g in gb.items():
        if g > G:
            continue
        b = bend(v)
        if b == 0:
            continue  # the two planes are drawn separately
        x, y, z, r = geom(v)
        assert abs(r - 1 / b) < 1e-9
        rows.append((g, gen_a(v), b, x, y, z, r))
    rows.sort()
    # sanity: strategy-(a) counts by generation, among spheres with g_b <= G
    from collections import Counter
    ca = Counter(r[1] for r in rows)
    cb = Counter(r[0] for r in rows)
    print("g_b counts", [cb[k] for k in range(G + 1)])
    print("g_a counts (complete up to G when g_b <= g_a)", [ca[k] for k in range(G + 1)])
    print("pairs with g_a < g_b:", sum(1 for r in rows if r[1] < r[0]))
    buf = b"".join(struct.pack("<ffffBBH", x, y, z, r, g, ga, min(b, 65535)) for g, ga, b, x, y, z, r in rows)
    json.dump({"G": G, "n": len(rows), "stride": 20, "data": base64.b64encode(buf).decode()},
              open(f"results/viz_spheres_{G}.json", "w"))
    print(len(rows), "spheres,", len(buf), "bytes")
