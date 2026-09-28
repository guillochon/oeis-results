\\ What does ellrankinit contain? Is the cubic field bnf inside it (so we can certify it)?
default(nbthreads, 1);
{
  my(v = 11244925, E = ellinit([0, 0, 0, 0, v]), R = ellrankinit(E));
  print("type ", type(R), " length ", #R);
  for (i = 1, #R, print(i, ": ", type(R[i]), if (type(R[i]) == "t_VEC", Str(" len ", #R[i]), "")));
  \\ search one level down for a bnf (a t_VEC of length 10 with a nf in position 7)
  for (i = 1, #R, if (type(R[i]) == "t_VEC",
    for (j = 1, #R[i], my(z = R[i][j]); if (type(z) == "t_VEC" && #z == 10, print("candidate bnf at [", i, "][", j, "]: pol ", z[7][1])))));
}
