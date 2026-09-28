\\ h3_cubic_test.gp -- unconditional 3-ranks from cubic fields (Hasse: a fundamental D has
\\ (3^r - 1)/2 cubic fields of discriminant D, r = 3-rank of Cl(Q(sqrt D))).
\\ Compare with quadclassunit (GRH) on a small range and time nflist("S3", range).
\\   gp -q selmer3_bound.gp h3_cubic_test.gp
\\ Result (2026-09-27): 0 disagreements for |D| <= 2*10^5, 17 s.
default(parisizemax, 2*10^9); default(threadsizemax, 2*10^9);
{
\\ count cubic fields per discriminant (both signs) and compare with quadclassunit
my(X = 2*10^5, cnt = Map(), bad = 0);
foreach(concat(nflist("S3", [1, X], 0), nflist("S3", [1, X], 1)), P,
  my(d = nfdisc(P)); if(isfundamental(d), mapput(cnt, d, if(mapisdefined(cnt, d), mapget(cnt, d), 0) + 1)));
for(d = -X, X, if(!isfundamental(d) || d == 1, next);
  my(n = if(mapisdefined(cnt, d), mapget(cnt, d), 0), r = logint(2*n + 1, 3));
  if(r != h3(d), bad++));
print("fundamental |D| <= ", X, ": disagreements between cubic-field count and quadclassunit: ", bad);
}
