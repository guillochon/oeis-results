\\ selmer3_bound.gp -- cheap upper bound on rank(y^2 = x^3 + k) via descent by 3-isogeny.
\\
\\ Source: Cohen & Pazuki, "Elementary 3-descent with a 3-isogeny", Acta Arith. 140 (2009),
\\ arXiv:0903.4963, with a = 0 (Mordell curves).
\\   Write k = D*b^2, D a fundamental discriminant (D = 1 allowed), K = Q(sqrt D).
\\   Prop. 2.2: |Im alpha| * |Im alpha'| = 3^(r + delta), delta = 1 iff D in {1, -3}.
\\   Cor. 4.3 + Lemma 4.4: Im alpha sits in a group of F_3-dimension at most
\\        h3(K) + dim U(K)/U(K)^3 + #{primes p | 2b split in K}
\\   (for D = 1, i.e. K = Q x Q, Thm 3.1(3): at most omega(2b)).
\\   The dual curve is y^2 = x^3 - 27k (same rank).
\\ So  rank <= s(k) + s(-27k) - delta.  This ignores all local conditions, so it is an
\\ upper bound on the 3-isogeny Selmer bound, which is itself an upper bound on the rank.
\\
\\ WARNING: quadclassunit is conditional on GRH. For an unconditional run, h3 must come from
\\ a certified source (cubic-field enumeration, see README). Here it is only used for timing
\\ and for testing the formula.

h3(D) = if(D == 1, 0, my(c = quadclassunit(D).cyc); sum(i = 1, #c, c[i] % 3 == 0));

\\ F_3-dimension of the ambient group for the 3-descent map on y^2 = x^3 + v
side(v) =
{
  my(D = coredisc(v), b2 = sqrtint(4*v/D), f = factor(b2)[,1], s);
  if(D == 1, return(#f));                         \\ K = Q x Q
  s = h3(D) + (D > 0 || D == -3);                 \\ class group + units mod cubes
  s + sum(i = 1, #f, kronecker(D, f[i]) == 1);    \\ split primes dividing 2b
}

B3(k) = my(D = coredisc(k)); side(k) + side(-27*k) - (D == 1 || D == -3);

\\ is k sixth-power free?
sfree6(k) = my(f = factor(k)[,2]); for(i = 1, #f, if(f[i] >= 6, return(0))); 1;
