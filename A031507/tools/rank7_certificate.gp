\\ Certificate that y^2 = x^3 + 47550317 and y^2 = x^3 - 56877643 have rank exactly 7.
\\ Lower bound: 7 rational points (checked on the curve) with nonzero regulator, computed at two
\\ precisions (the height pairing matrix is positive definite iff the points are independent).
\\ Upper bound: the unconditional 3-isogeny bound B2 = min(2 sa - g, 2 sb + g) - d = 7, from
\\ the certified h3 table (sa, sb from the Rust filter) and Cassels' formula (g).
default(nbthreads, 1);
read("tools/stage2.gp");
{
  foreach ([[47550317, 1], [56877643, -1]], c,
    my(k = c[1], v = c[2] * k, E = ellinit([0, 0, 0, 0, v]), z = ellrank(E, 3), P = z[4]);
    printf("v = %d: ellrank %s, %d points\n", v, z[1..3], #P);
    for (i = 1, #P, if (!ellisoncurve(E, P[i]), error("point not on curve")));
    foreach ([38, 115], pr,
      default(realprecision, pr);
      my(M = ellheightmatrix(E, P), d = matdet(M));
      printf("   precision %d: regulator %.12g, min eigenvalue %.6g\n", pr, d, vecmin(qfjacobi(M)[1])));
    printf("   Cassels g = %d\n", cassels_g(v));
    print("   points: ", P));
}
