\\ calibrate.gp -- test and time selmer3_bound.gp.  Run from WSL:
\\   gp -q selmer3_bound.gp calibrate.gp
default(parisizemax, 2*10^9);

\\ 1. Record curves: the bound must be >= the rank.  (k, rank); negative k means y^2 = x^3 - |k|.
{rec = [[1,0],[2,1],[15,2],[113,3],[2089,4],[66265,5],[1358556,6],[47550317,7],[1632201497,8],
       [-1,0],[-2,1],[-11,2],[-174,3],[-2351,4],[-28279,5],[-975379,6],[-56877643,7],[-2520963512,8]];}
print("record curves: k, rank, B3");
for(i = 1, #rec, my(k = rec[i][1]); print(k, "  ", rec[i][2], "  ", B3(k), if(B3(k) < rec[i][2], "   *** FAIL", "")));

\\ 2. Small k: bound vs. the proven lower bound from ellrank (points found).
bad = 0; slack = vector(12);
{
for(k = -3000, 3000, if(k == 0 || !sfree6(abs(k)), next);
  my(r = ellrank(ellinit([0, k]))[1], b = B3(k));
  if(b < r, bad++; print("FAIL k=", k, " r>=", r, " B3=", b));
  slack[min(b - r, 11) + 1]++);
}
print("small k: failures = ", bad, "; histogram of B3 - (ellrank lower bound), 0..11: ", slack);

\\ 3. Timing and distribution near the top of the range (10^4 values each sign).
{
foreach([1, -1], sg,
  my(h = vector(15), t = getabstime(), n = 0);
  for(k = 47540000, 47550000, if(!sfree6(k), next); n++; h[min(B3(sg*k), 14) + 1]++);
  print("sign ", sg, ": ", n, " curves in ", getabstime() - t, " ms; histogram of B3 (0..14): ", h));
}
