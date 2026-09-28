\\ Stage 4: close survivors with witnesses (unconditional).
\\ Input: infile (lines "k sa sb d g"), sgn, outfile, PMAX (local tests at primes of S up to PMAX).
\\ Output lines "k A B ma mb newbound": A = 2sa-g-d, B = 2sb+g-d are the stage-2 terms; ma, mb
\\ are certified witness counts for side a (u = -27v, field Q(sqrt(-3v))) and side b (u = v).
\\ newbound = min(2(sa-ma)-g, 2(sb-mb)+g) - d, a proven upper bound on the rank.
default(nbthreads, 1);
read("tools/witness.gp");
Sof(k) = setunion([2, 3], Set(factor(k)[, 1]~));
{
  my(L = readstr(infile), rows = vector(#L), hist = vector(12));
  for (i = 1, #L,
    my(w = apply(eval, strsplit(L[i], " ")), k = w[1], sa = w[2], sb = w[3], d = w[4], g = w[5]);
    my(v = sgn * k, S = Sof(k), P = select(p -> p <= PMAX, S));
    my(A = 2*sa - g - d, B = 2*sb + g - d, ma = 0, mb = 0, nb);
    \\ the smaller term decides; improve it, or either one on a tie (b first: smaller field)
    my(needb = (B - TGT + 1) \ 2, needa = (A - TGT + 1) \ 2);
    if (B <= A && !issquare(v), mb = witness_m(v, S, P, min(needb, 2))[1]);
    nb = min(2*(sa - ma) - g, 2*(sb - mb) + g) - d;
    if (nb > TGT && A <= B && !issquare(-3*v), ma = witness_m(-27*v, S, P, min(needa, 2))[1]);
    nb = min(2*(sa - ma) - g, 2*(sb - mb) + g) - d;
    hist[nb + 1]++;
    rows[i] = Str(k, " ", A, " ", B, " ", ma, " ", mb, " ", nb));
  write(outfile, strjoin(rows, "\n"));
  print(infile, ": ", #L, " curves; new bound histogram ", hist);
}
