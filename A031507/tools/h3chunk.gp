\\ Unconditional 3-ranks of quadratic class groups from cubic-field counts (Hasse):
\\ for a fundamental discriminant D, #{cubic fields of discriminant D} = (3^h3(D) - 1)/2.
\\ nflist("S3", [lo, hi]) enumerates (unconditionally, Belabas' method) all non-Galois cubic
\\ fields with lo <= |disc| <= hi, both signs; cyclic cubic fields have square discriminant,
\\ never fundamental. Input variables: lo, hi, out (file name).
\\ Writes lines "D h3" for every fundamental D with h3 >= 1, and a summary line to out.done.
default(nbthreads, 1);
{
  my(L = nflist("S3", [lo, hi]), v = List(), bad = 0, cnt = vector(8));
  for (i = 1, #L, my(d = nfdisc(L[i])); if (isfundamental(d), listput(v, d)));
  v = vecsort(Vec(v));
  my(i = 1, rows = List());
  while (i <= #v,
    my(j = i); while (j < #v && v[j + 1] == v[i], j++);
    my(n = j - i + 1, t = 2*n + 1, h = valuation(t, 3));
    if (3^h != t, bad++; print("NOT A POWER OF 3: D=", v[i], " fields=", n));
    cnt[h] += 1;
    listput(rows, Str(v[i], " ", h));
    i = j + 1);
  \\ one write per chunk: per-line writes to /mnt/c are very slow under WSL
  if (#rows, write(out, strjoin(Vec(rows), "\n")));
  write(Str(out, ".done"), "lo=", lo, " hi=", hi, " fields=", #L, " fundamental_fields=", #v, " D_by_h3=", cnt, " bad=", bad);
}
