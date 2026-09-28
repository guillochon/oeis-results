\\ Stage 6: unconditional 2-descent for the curves the 3-isogeny method could not close.
\\ ellrankinit(E) holds the bnf of the cubic field Q[x]/(x^3 + v) used by the 2-descent
\\ (component [3][1]). bnfcertify proves that bnf (class group and units) without GRH; ellrank on
\\ the same structure then gives an upper bound R that no longer depends on GRH.
\\ Input: infile (one k per line), sgn, outfile. Output: "k certified R r ms".
default(nbthreads, 1);
{
  my(L = readstr(infile), rows = vector(#L), bad = 0, maxR = 0);
  for (i = 1, #L,
    my(k = eval(L[i]), v = sgn * k, t = getabstime());
    my(E = ellinit([0, 0, 0, 0, v]), RI = ellrankinit(E), bnf = RI[3][1]);
    if (bnf.pol != x^3 + v && poldegree(bnf.pol) != 3, error("unexpected field for v = ", v));
    my(ok = bnfcertify(bnf), z = ellrank(RI));
    if (!ok, bad++);
    if (ok, maxR = max(maxR, z[2]));
    rows[i] = Str(k, " ", ok, " ", z[2], " ", z[1], " ", getabstime() - t));
  write(outfile, strjoin(rows, "\n"));
  print(infile, ": ", #L, " curves; certification failures ", bad, "; max certified R ", maxR);
}
