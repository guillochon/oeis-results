\\ The a(6) curves have rank exactly 6: ellrank points (lower bound) and the certified
\\ 2-descent (upper bound, bnfcertify on the cubic field used by ellrank).
default(nbthreads, 1);
{
  foreach ([1358556, -975379], v,
    my(E = ellinit([0, 0, 0, 0, v]), RI = ellrankinit(E), ok = bnfcertify(RI[3][1]), z = ellrank(RI, 3));
    default(realprecision, 60);
    printf("v = %d: certified %d, ellrank [r, R] = [%d, %d], regulator %.10g\n", v, ok, z[1], z[2], matdet(ellheightmatrix(E, z[4]))));
}
