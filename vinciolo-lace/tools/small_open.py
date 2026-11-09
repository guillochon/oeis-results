"""Open 8-connected regions by size: panels vs density-matched random no-pinhole designs."""
import numpy as np, collections
from growth_rules import P, N, comp_sizes, S8
tot = 1.0; logp = 0
for k, A in P.items():
    s = sorted(comp_sizes(A, 0, S8)); small = [x for x in s if x < 12]
    nb = [sorted(x for x in comp_sizes(B, 0, S8) if x < 5) for B in N[k]]
    rate = np.mean([len(x) == 0 for x in nb]); mean_small = np.mean([len(x) for x in nb])
    tot *= max(rate, 1/40)
    print(f'{k}  open regions {len(s):3d}  smallest {s[:6]}  random: no region<5 in {rate:5.0%}, mean #<5 = {mean_small:.1f}')
print('product of per-panel random pass rates:', tot)
c = collections.Counter(x for k in P for B in N[k] for x in comp_sizes(B, 0, S8) if x < 12)
print('random size histogram (<12):', sorted(c.items()))
c = collections.Counter(x for A in P.values() for x in comp_sizes(A, 0, S8) if x < 30)
print('panel size histogram (<30):', sorted(c.items()))
