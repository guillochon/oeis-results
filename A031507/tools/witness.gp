\\ Witnesses that shrink a 3-isogeny Selmer group below its ambient bound (unconditional).
\\
\\ Setting: E_u: y^2 = x^3 + u, K = Q(sqrt u) (u not a square), S = primes dividing 6k.
\\ The descent map E_u(Q) -> K*/K*^3, (x, y) -> y + sqrt(u), lands in the Selmer group of the
\\ 3-isogeny onto E_u, inside the minus part A of K(S,3) (dim A <= s). For c = c1 + c2*sqrt(u)
\\ with N(c) a rational cube, the class of c is in the Selmer group iff the covering
\\     C_c : Z^3 = F_c(X, Y) = c2 X^3 + 3 c1 X^2 Y + 3 u c2 X Y^2 + u c1 Y^3
\\ (from c (X + Y sqrt u)^3 - conj(c) (X - Y sqrt u)^3 = 2 sqrt(u) Z^3) has Q_p-points for all p.
\\ If c_1..c_m are explicit elements of A and every nontrivial product c_1^e1...c_m^em is locally
\\ insoluble at some prime, then W = <c_i> meets Sel trivially, so dim Sel <= dim A - m <= s - m.
\\ (A trivial product would be soluble, so this also certifies independence.)
\\ bnfinit (GRH) is used only to FIND candidate elements; every element is verified exactly.

\\ ---------- local solubility of Z^3 = F(X,Y) over Q_p ----------
\\ is the p-adic unit a (an integer prime to p) a cube in Z_p?
iscubeunit(a, p) =
{
  if (p == 3, my(r = a % 9); return(r == 1 || r == 8));
  if (p % 3 == 2, return(1));
  Mod(a, p)^((p - 1) / 3) == 1;
}

\\ does g(t) (integer polynomial in t) take a p-adic cube value (or 0) for some t in Z_p?
\\ returns 1 (yes), 0 (no), -1 (undecided: depth exhausted)
locsol_poly(g, p, depth) =
{
  my(c = Vecrev(g), v0, m = oo, need = if (p == 3, 2, 1));
  if (g == 0, return(1));
  if (c[1] == 0, return(1));
  v0 = valuation(c[1], p);
  for (i = 2, #c, if (c[i] != 0, m = min(m, valuation(c[i], p))));
  if (m >= v0 + need,
    return(v0 % 3 == 0 && iscubeunit(c[1] / p^v0, p)));
  if (depth <= 0, return(-1));
  my(res = 0);
  for (j = 0, p - 1,
    my(r = locsol_poly(subst(g, 't, j + p * 't), p, depth - 1));
    if (r == 1, return(1));
    if (r == -1, res = -1));
  res;
}

\\ Fast version for p != 3: descend only into residues that are roots of the reduction.
\\ g = p^a h with h primitive. For t with h(t) != 0 mod p: v(g(t)) = a and the cube class of the
\\ unit part is that of h(t) mod p. Such t give a cube iff a = 0 mod 3 and (p = 2 mod 3, or h(t)
\\ is a nonzero cube mod p). Returning 1 ("soluble") is always safe for our use; a 0 is returned
\\ only after every residue class has been ruled out exactly.
locsol_fast(g, p, depth) =
{
  if (g == 0, return(1));
  my(a = valuation(content(g), p), h = g / p^a, hb = Mod(h, p), R);
  if (poldegree(lift(hb)) <= 0, R = [], R = polrootsmod(lift(hb), p));
  \\ residues t mod p that are not roots of hb
  if (a % 3 == 0 && p - #R > 0,
    if (p % 3 == 2, return(1));
    \\ p = 1 mod 3: is some non-root value a cube?
    if (p <= 50,
      for (t = 0, p - 1, my(x = subst(lift(hb), 't, t) % p); if (x && Mod(x, p)^((p-1)/3) == 1, return(1))),
      \\ p > 50: unless hb = const * (linear)^3 or const, Weil's bound gives a cube value
      my(hh = lift(hb), dg = poldegree(hh), cst = 0);
      if (dg <= 0, cst = hh,
        my(f = factormod(hh, p));
        if (#f~ == 1 && f[1, 2] == 3 && poldegree(f[1, 1]) == 1, cst = pollead(hh)));
      if (cst == 0, return(1));
      if (Mod(cst, p)^((p-1)/3) == 1, return(1))));
  if (depth <= 0, return(if (#R, -1, 0)));
  my(res = 0);
  for (i = 1, #R,
    my(j = lift(R[i]), r = locsol_fast(subst(g, 't, j + p * 't), p, depth - 1));
    if (r == 1, return(1));
    if (r == -1, res = -1));
  res;
}

\\ 1 if Z^3 = F(X,Y) has a Q_p-point, 0 if provably none, -1 if undecided.
\\ Charts: (X, 1) with X in Z_p, and (1, pY) with Y in Z_p; together they cover P^1(Q_p).
locsol(F, p) =
{
  my(fx = subst(subst(F, 'Y, 1), 'X, 't), fy = subst(subst(F, 'X, 1), 'Y, p * 't));
  my(sol(g) = if (p == 3, locsol_poly(g, 3, 40), locsol_fast(g, p, 60)));
  my(r1 = sol(fx));
  if (r1 == 1, return(1));
  my(r2 = sol(fy));
  if (r2 == 1, return(1));
  if (r1 == -1 || r2 == -1, return(-1));
  0;
}

\\ covering form for c = c1 + c2 sqrt(u), c1, c2 rational; scaled by a rational cube to integers
covering(c1, c2, u) =
{
  my(F = c2*'X^3 + 3*c1*'X^2*'Y + 3*u*c2*'X*'Y^2 + u*c1*'Y^3, den = denominator(content(F)));
  my(l = 1, fd = factor(den));
  for (i = 1, #fd~, l *= fd[i, 1]^ceil(fd[i, 2] / 3));
  F * l^3;
}

\\ first prime of P at which the covering of c = c1 + c2 sqrt u is provably insoluble, else 0
insoluble_at(c1, c2, u, P) =
{
  my(F = covering(c1, c2, u));
  for (i = 1, #P, if (locsol(F, P[i]) == 0, return(P[i])));
  0;
}

\\ ---------- candidate elements of A, each verified ----------
\\ K = Q(w), w^2 = D (fundamental), sqrt(u) = r*w. Element gam (polmod in w) -> c = gam/conj(gam),
\\ N(c) = 1. Returns [c1, c2, why] (c = c1 + c2 sqrt u) if c is verified to lie in K(S,3), else 0.
addc(bnf, D, r, S, gam, why) =
{
  my(ga = Mod(lift(gam), 'w^2 - D), cg = Mod(subst(lift(ga), 'w, -'w), 'w^2 - D), c = ga / cg);
  if (norm(c) != 1, error("bad norm ", why));
  my(f = factor(abs(norm(ga))));
  for (i = 1, #f~,
    my(q = f[i, 1]);
    if (setsearch(S, q), next);
    my(pr = idealprimedec(bnf, q));
    for (j = 1, #pr, if (idealval(bnf, lift(c), pr[j]) % 3 != 0, return(0))));
  my(z = lift(c));
  [polcoef(z, 0, 'w), polcoef(z, 1, 'w) / r, why];
}

\\ list of verified candidates [c1, c2, why] in the minus part of K(S,3), K = Q(sqrt u)
cands(u, S) =
{
  my(D = coredisc(u), r2 = u / D, r, bnf, L = List(), z);
  if (!issquare(r2, &r), error("u/D not a square for u = ", u));
  bnf = bnfinit('w^2 - D, 1);
  if (D > 0, z = addc(bnf, D, r, S, Mod(nfbasistoalg(bnf, bnf.fu[1]), 'w^2 - D), "unit"); if (z, listput(L, z)));
  foreach (S, p,
    if (kronecker(D, p) != 1, next);
    my(pr = idealprimedec(bnf, p)[1], e = bnfisprincipal(bnf, pr, 0), o = 1);
    for (i = 1, #e, if (e[i], o = lcm(o, bnf.cyc[i] / gcd(bnf.cyc[i], e[i]))));
    my(I = idealpow(bnf, pr, o), zz = bnfisprincipal(bnf, I, 3));
    if (zz[1] != 0, next);
    my(gam = nfbasistoalg(bnf, zz[2]));
    if (idealhnf(bnf, gam) != idealhnf(bnf, I), error("principal check failed"));
    z = addc(bnf, D, r, S, gam, Str("split ", p)); if (z, listput(L, z)));
  for (i = 1, #bnf.cyc,
    if (bnf.cyc[i] % 3, next);
    my(a = idealpow(bnf, bnf.gen[i], bnf.cyc[i] / 3), I = idealpow(bnf, a, 3), zz = bnfisprincipal(bnf, I, 3));
    if (zz[1] != 0, error("a^3 not principal?"));
    my(gam = nfbasistoalg(bnf, zz[2]));
    if (idealhnf(bnf, gam) != idealhnf(bnf, I), error("principal check failed"));
    z = addc(bnf, D, r, S, gam, Str("class ", i)); if (z, listput(L, z)));
  Vec(L);
}

\\ multiply elements (c1 + c2 sqrt u)
cmul(a, b, u) = [a[1] * b[1] + u * a[2] * b[2], a[1] * b[2] + a[2] * b[1]];
cpow(a, e, u) = my(r = [1, 0]); for (i = 1, e, r = cmul(r, a, u)); r;

\\ Full-span search. Enumerates every nonzero combination prod C_i^e_i (e up to scalar), tests
\\ it at all primes of S, and returns the largest certified m <= mmax (<= 2) with an m-dim
\\ subspace all of whose nonzero elements are locally insoluble. Returns [m, #span, #insoluble].
witness_span(u, S, mmax) =
{
  my(C = cands(u, S), n = #C, bad = Map(), nb = 0, tot = 0);
  if (n == 0, return([0, 0, 0]));
  my(pw = vector(n, i, [[1, 0], C[i], cmul(C[i], C[i], u)]));
  forvec (e = vector(n, i, [0, 2]),
    my(f = 0); for (i = 1, n, if (e[i], f = e[i]; break));
    if (f != 1, next);                                         \\ projective representative
    my(z = [1, 0]); for (i = 1, n, if (e[i], z = cmul(z, pw[i][e[i] + 1], u)));
    tot++;
    if (insoluble_at(z[1], z[2], u, S), mapput(bad, e, 1); nb++));
  if (nb == 0, return([0, tot, 0]));
  if (mmax <= 1, return([1, tot, nb]));
  my(B = Vec(Mat(bad)[, 1]));
  my(norm3(e) = my(f = 0); for (i = 1, #e, if (e[i] % 3, f = e[i] % 3; break)); if (f == 2, e = -e); e % 3);
  for (i = 1, #B, for (j = i + 1, #B,
    my(x = B[i], y = B[j]);
    if (mapisdefined(bad, norm3(x + y)) && mapisdefined(bad, norm3(x + 2 * y)), return([2, tot, nb]))));
  [1, tot, nb];
}

\\ Largest certified m <= mmax (<= 2): m elements whose span meets Sel trivially, using local
\\ tests at the primes P. Returns [m, descriptions].
witness_m(u, S, P, mmax) =
{
  my(C = cands(u, S), bad = List());
  for (i = 1, #C, if (insoluble_at(C[i][1], C[i][2], u, P), listput(bad, C[i])));
  \\ products c_i c_j^e can be insoluble even when the single elements are soluble
  if (#bad == 0,
    for (i = 1, #C, for (j = i + 1, #C, for (e = 1, 2,
      my(z = cmul(C[i], cpow(C[j], e, u), u));
      if (insoluble_at(z[1], z[2], u, P), listput(bad, [z[1], z[2], Str(C[i][3], "*", C[j][3], "^", e)]))))));
  if (#bad == 0, return([0, []]));
  if (mmax == 1, return([1, [bad[1][3]]]));
  \\ m = 2: need x, y with x, y, x*y, x*y^2 all insoluble (x^2 etc. are inverses)
  my(B = Vec(bad));
  for (i = 1, #B, for (j = i + 1, #B,
    my(x = B[i], y = B[j], ok = 1);
    for (e = 1, 2, my(z = cmul(x, cpow(y, e, u), u)); if (!insoluble_at(z[1], z[2], u, P), ok = 0; break));
    if (ok, return([2, [x[3], y[3]]]))));
  [1, [B[1][3]]];
}
