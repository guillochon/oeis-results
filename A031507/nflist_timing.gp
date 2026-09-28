\\ nflist_timing.gp -- how fast does PARI enumerate cubic (S3) fields by discriminant range?
\\   gp -q nflist_timing.gp
default(threadsizemax, 2*10^9); default(parisizemax, 4*10^9);
{
foreach([10^5, 10^6, 4*10^6], X,
  my(t = getabstime(), w = getwalltime(), L = nflist("S3", [1, X]));
  print("|disc| <= ", X, ": ", #L, " cubic fields; cpu ", getabstime() - t, " ms, wall ", getwalltime() - w, " ms"));
}
