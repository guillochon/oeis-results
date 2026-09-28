\\ Stage 2: refine the 3-isogeny bound with Cassels' formula (unconditional).
\\
\\ For E = E_v: y^2 = x^3 + v and E' = E_{-27v}, with phi: E -> E' of degree 3:
\\   |Sel^phi(E)| / |Sel^phihat(E')| = |E(Q)[phi]| Omega(E') prod c_p(E') / (|E'(Q)[phihat]| Omega(E) prod c_p(E))
\\ (Cassels 1965). With ellbsd(E) = Omega(E) prod c_p(E) / |E(Q)_tors|^2 (minimal model, Tate's
\\ algorithm; no L-functions) this gives g = dim Sel^phi - dim Sel^phihat exactly.
\\ With dim Sel^phi <= sa (ambient bound in Q(sqrt(-3v))) and dim Sel^phihat <= sb (in Q(sqrt v)):
\\   rank <= dim Sel^phi + dim Sel^phihat - delta <= min(2*sa - g, 2*sb + g) - delta.

\\ exact g (and a check that the Cassels ratio is a power of 3)
cassels_g(v) =
{
  my(E = ellinit([0, 0, 0, 0, v]), F = ellinit([0, 0, 0, 0, -27*v]));
  my(tE = elltors(E)[1], tF = elltors(F)[1]);
  my(kE = if (issquare(v), 3, 1), kF = if (issquare(-27*v), 3, 1));
  my(q = kE * ellbsd(F) * tF^2 / (kF * ellbsd(E) * tE^2));
  my(g = round(log(q) / log(3)));
  if (abs(q / 3.^g - 1) > 1e-20, error("Cassels ratio not a power of 3 for v = ", v, ": ", q));
  g;
}

B2(sa, sb, d, g) = min(2*sa - g, 2*sb + g) - d;
