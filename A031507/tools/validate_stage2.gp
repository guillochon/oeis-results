\\ Validation of the stage-2 bound: for every sixth-power-free k <= K, both signs, check that
\\ B2 >= r (the proven lower bound on the rank from ellrank's independent points), and count
\\ how often the bound with the Cassels ratio inverted would fail (a control).
\\ Input: data/cand_{plus,minus}_0_K.txt from `a031507 filter K 0`.  Usage: K = ...; read(...)
default(nbthreads, 1);
default(realprecision, 60);
read("tools/stage2.gp");
{
  my(bad = 0, badinv = 0, n = 0, hist = vector(12), histB = vector(12));
  foreach ([["plus", 1], ["minus", -1]], P,
    my(L = readstr(Str("data/cand_", P[1], "_0_", K, ".txt")));
    for (i = 1, #L,
      my(w = apply(eval, strsplit(L[i], " ")), k = w[1], sa = w[2], sb = w[3], d = w[4], v = P[2]*k);
      my(g = cassels_g(v), b2 = B2(sa, sb, d, g), binv = B2(sa, sb, d, -g));
      my(r = ellrank(ellinit([0, 0, 0, 0, v]))[1]);
      n++; hist[b2 + 1]++; histB[sa + sb - d + 1]++;
      if (r > b2, bad++; print("VIOLATION v=", v, " r=", r, " B2=", b2, " sa=", sa, " sb=", sb, " g=", g));
      if (r > binv, badinv++)));
  print("curves ", n, "; B2 < r: ", bad, "; inverted-ratio bound < r: ", badinv);
  print("B histogram  ", histB);
  print("B2 histogram ", hist);
}
