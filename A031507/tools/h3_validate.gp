\\ Compare the unconditional table (data/h3, |D| <= M) with quadclassunit (GRH) for every
\\ fundamental D with 1 < |D| <= M.  Usage: M = ...; read("tools/h3_validate.gp")
default(nbthreads, 1);
T = Map();
{
  my(fs = externstr("ls data/h3/h3_*.txt"));
  for (i = 1, #fs,
    my(L = readstr(fs[i]));
    for (j = 1, #L, my(w = strsplit(L[j], " ")); mapput(T, eval(w[1]), eval(w[2]))));
  my(mis = 0, tot = 0, h = 0);
  forstep (D = -M, M, 1,
    if (D != 1 && D != 0 && isfundamental(D),
      tot++;
      my(t = 0, c = quadclassunit(D).cyc);
      mapisdefined(T, D, &t);
      h = sum(i = 1, #c, c[i] % 3 == 0);
      if (h != t, mis++; if (mis < 10, print("mismatch D=", D, " table=", t, " qcu=", h)))));
  print("validated ", tot, " fundamental D with |D| <= ", M, ": mismatches ", mis);
}
