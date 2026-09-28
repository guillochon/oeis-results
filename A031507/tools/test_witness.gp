\\ Tests for witness.gp.
\\ (1) soundness of the local test: images y + sqrt(u) of rational points must be soluble at every
\\     prime of S (they have a rational point on the covering).
\\ (2) the improved bound never drops below the proven lower bound r of ellrank, for all
\\     sixth-power-free 1 <= k <= K, both signs.
default(nbthreads, 1);
default(realprecision, 60);
read("tools/stage2.gp");
read("tools/witness.gp");
Sof(k) = setunion([2, 3], Set(factor(k)[, 1]~));
smallP(S) = select(p -> p <= 200, S);
{
  my(bad1 = 0, tested = 0);
  for (k = 2, 400, foreach ([1, -1], sg,
    my(u = sg * k);
    if (issquare(u) || !issquarefree(k), next);
    my(E = ellinit([0, 0, 0, 0, u]), L = ellrank(E)[4], S = Sof(k));
    for (i = 1, #L, my(P = L[i]);
      tested++;
      my(D = coredisc(u), r2 = u / D);
      \\ c = y + sqrt(u): c1 = y, c2 = 1; N(c) = y^2 - u = x^3, a rational cube
      if (insoluble_at(P[2], 1, u, smallP(S)), bad1++; print("UNSOUND: u=", u, " P=", P)))));
  print("test 1: ", tested, " point images tested, false insolubility: ", bad1);
}
{
  my(viol = 0, n = 0, gain = vector(10));
  for (k = 1, K, foreach ([1, -1], sg,
    my(v = sg * k, f = factor(k), ok = 1);
    for (i = 1, #f~, if (f[i, 2] >= 6, ok = 0));
    if (!ok, next);
    my(S = Sof(k), P = smallP(S));
    \\ ambient bounds from the unconditional formula, with GRH-free h3 replaced here by quadclassunit
    \\ (only for this test; the production run reads the certified table)
    my(h3q(D) = if (D == 1, 0, my(c = quadclassunit(D).cyc); sum(i = 1, #c, c[i] % 3 == 0)));
    my(sfield(D) = if (D == 1, #S, h3q(D) + (D > 0 || D == -3) + sum(i = 1, #S, kronecker(D, S[i]) == 1)));
    my(Da = coredisc(-3 * v), Db = coredisc(v), sa = sfield(Da), sb = sfield(Db), d = (Da == 1 || Db == 1));
    my(g = cassels_g(v), A = 2 * sa - g - d, B = 2 * sb + g - d, ma = 0, mb = 0);
    if (Db != 1, mb = witness_m(v, S, P, 2)[1]);
    if (Da != 1, ma = witness_m(-27 * v, S, P, 2)[1]);
    my(bnew = min(2 * (sa - ma) - g, 2 * (sb - mb) + g) - d, r = ellrank(ellinit([0, 0, 0, 0, v]))[1]);
    n++; gain[min(A, B) - bnew + 1]++;
    if (r > bnew, viol++; print("VIOLATION v=", v, " r=", r, " new bound=", bnew, " sa sb g ma mb=", [sa, sb, g, ma, mb]))));
  print("test 2: ", n, " curves; bound below r: ", viol, "; reduction histogram (0,1,2,3..) ", gain);
}
