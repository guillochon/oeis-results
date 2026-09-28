\\ How long does bnfcertify take on the 2-descent cubic field of a leftover curve?
default(nbthreads, 1);
{
  foreach ([5814352, 11244925, -3170339], v,
    my(t = getabstime(), bnf = bnfinit(x^3 + v, 1), t1 = getabstime() - t);
    printf("v=%d: disc %s, bnfinit %d ms, cyc %s\n", v, factor(nfdisc(x^3 + v)), t1, bnf.cyc);
    t = getabstime();
    my(ok = bnfcertify(bnf));
    printf("   bnfcertify = %d in %d ms\n", ok, getabstime() - t));
}
