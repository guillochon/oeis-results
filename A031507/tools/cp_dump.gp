\\ Independent cross-check: the agent's Cohen-Pazuki bound B3 (h3 from quadclassunit, GRH),
\\ for 1 <= k <= K, both signs. Output lines "k B3(+k) B3(-k)", -1 if not sixth-power free.
\\ Usage: K = ...; read("tools/cp_dump.gp")
default(nbthreads, 1);
read("selmer3_bound.gp");
{
  my(rows = vector(K));
  for (k = 1, K,
    if (!sfree6(k), rows[k] = Str(k, " -1 -1"),
      rows[k] = Str(k, " ", B3(k), " ", B3(-k))));
  write(Str("data/cpdump_", K, ".txt"), strjoin(rows, "\n"));
}
