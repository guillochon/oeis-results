\\ Check the PARI lines of oeis_submission.md against the b-files (n = 0..1000).
\\ Run from attacks/A154640:  gp -q check_pari.gp
default(parisizemax, 2000000000);
a639(n) = polcoef((1+2*x-x^2+2*x^3)/(1-3*x-6*x^2+12*x^3) + O(x^(n+1)), n);
a640(n) = polcoef(5*(1-x)*(1-x-3*x^2)/(1-3*x-3*x^2+6*x^3) + O(x^(n+1)), n);
a641(n) = polcoef((1-x)^2*(1+x)*(1-x-3*x^2)/(1-5*x-8*x^2+21*x^3+9*x^4-6*x^5) + O(x^(n+1)), n);
a642(n) = polcoef((1-x)^3*(1+x)*(1-x-3*x^2)/(1-10*x-7*x^2+85*x^3-24*x^4-27*x^5+6*x^6) + O(x^(n+1)), n);
chk(f, file) = {
  my(L = readstr(file), bad = 0, s);
  for(i = 1, #L, s = strsplit(L[i], " "); if(f(eval(s[1])) != eval(s[2]), bad++));
  printf("%s: %d terms, %d mismatches\n", file, #L, bad);
}
chk(a639, "bfiles/b154639.txt");
chk(a640, "bfiles/b154640.txt");
chk(a641, "bfiles/b154641.txt");
chk(a642, "bfiles/b154642.txt");
quit;
