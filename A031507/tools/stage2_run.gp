\\ Apply the stage-2 (Cassels) bound to a candidate file.
\\ Input: infile (lines "k sa sb d"), sign (+1/-1), outfile.
\\ Writes "k sa sb d g B2" for every line to outfile, and prints a B2 histogram.
default(nbthreads, 1);
default(realprecision, 60);
read("tools/stage2.gp");
{
  my(L = readstr(infile), rows = vector(#L), hist = vector(16));
  for (i = 1, #L,
    my(w = apply(eval, strsplit(L[i], " ")), k = w[1], sa = w[2], sb = w[3], d = w[4]);
    my(g = cassels_g(sgn * k), b2 = B2(sa, sb, d, g));
    hist[b2 + 1]++;
    rows[i] = Str(k, " ", sa, " ", sb, " ", d, " ", g, " ", b2));
  if (#rows, write(outfile, strjoin(rows, "\n")));
  print(infile, ": ", #L, " curves; B2 histogram ", hist);
}
