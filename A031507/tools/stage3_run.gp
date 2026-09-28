\\ Stage 3 on the stage-2 output: parity, then ellrank.
\\ Input: infile (lines "k sa sb d g B2"), sgn, outfile.
\\ For B2 >= 7 only:
\\  - w = ellrootno(E) (unconditional). If d = 0 (no rational 3-torsion on E or E'), then
\\    rank <= s_inf <= B2 and (-1)^s_inf = w (3-parity theorem, Dokchitser-Dokchitser), so
\\    B3p = largest integer <= B2 with the parity of w. If d = 1 we keep B3p = B2.
\\  - for B3p >= 7: ellrank(E) -> [r, R] (R assumes GRH via bnfinit of the cubic field).
\\ Writes "k B2 w B3p r R" for every line with B2 >= 7.
default(nbthreads, 1);
{
  my(L = readstr(infile), rows = List(), nk = 0, nr = 0, n7 = 0, nR7 = 0);
  for (i = 1, #L,
    my(w6 = apply(eval, strsplit(L[i], " ")), k = w6[1], d = w6[4], b2 = w6[6]);
    if (b2 < 7, next);
    my(v = sgn * k, E = ellinit([0, 0, 0, 0, v]), w = ellrootno(E));
    my(b3 = if (d == 0 && (b2 % 2 == 0) != (w == 1), b2 - 1, b2));
    my(r = -1, R = -1);
    if (b3 >= 7, n7++; my(z = ellrank(E)); r = z[1]; R = z[2]; if (R >= 7, nR7++));
    listput(rows, Str(k, " ", b2, " ", w, " ", b3, " ", r, " ", R)));
  if (#rows, write(outfile, strjoin(Vec(rows), "\n")));
  print(infile, ": B2>=7: ", #rows, "; after parity B>=7: ", n7, "; ellrank R>=7: ", nR7);
}
