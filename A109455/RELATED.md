# Other OEIS sequences covered by the A109455 work

The n = 9 identity-class run of the A109455 engine enumerates every threshold function of 9 variables (one certified representative per NP class). Its by-products are listed below; the arithmetic ones are re-checked by `related_terms.py`. Logs: `logs/n9_identity.log`.

| Sequence | Result |
|---|---|
| [A000617](https://oeis.org/A000617) | **a(9) corrected** to 993061484 (the OEIS has 989913346); our enumeration, Kurz (2012) and A000619(9) agree |
| [A000609](https://oeis.org/A000609) | a(9) = 144130531453121108 **independently confirmed** (previously single-source) |
| [A002078](https://oeis.org/A002078) | **new a(9)** = 281814234754247 (counted directly; also the binomial mean of A000609) |
| [A001529](https://oeis.org/A001529) | **new a(9)** = 496618456, via a(n) = (A000617(n) + A001532(n))/2 |
| [A001530](https://oeis.org/A001530) | **new a(9)** = 495252138 = A001529(9) - A001529(8) |
| [A000615](https://oeis.org/A000615) | **new a(9)** = 143972777957019648 (inverse binomial transform of A000609) |
| [A002079](https://oeis.org/A002079) | **new a(9)** = 281196831947304 = A000615(9)/2^9 |
| [A001532](https://oeis.org/A001532) | a(9) = 175428 independently confirmed |
