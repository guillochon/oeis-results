\\ candidates_test.gp -- for k near the top of the range with B3(k) >= 7, time ellrank
\\ (2-descent; R is an unconditional upper bound on the rank).
\\   gp -q selmer3_bound.gp candidates_test.gp
default(parisizemax, 2*10^9);
{
foreach([1, -1], sg,
  for(k = 47500000, 47550000, if(!sfree6(k), next); my(b = B3(sg*k)); if(b < 7, next);
    my(t = getabstime(), R = ellrank(ellinit([0, sg*k])));
    print(sg*k, "  B3=", b, "  ellrank [r,R,s]=", R[1..3], "  ", getabstime() - t, " ms")));
}
