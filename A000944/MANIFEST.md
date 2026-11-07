# A000944 run manifest

Everything needed to replay the a(50) computation from this folder without network access, with
SHA-256 hashes. Regenerate with `make_manifest.py` in the source repository.

## Software

Python 3.12.0, numpy 2.2.6, scipy 1.16.0 (Windows 11, i5-12600K). networkx and nauty only for the brute-force checks.

## Run parameters (a(50))

* Box v <= 50, f <= 96 (`combine_full.py 50 96 10`).
* Sensed counts S: `sensed_p.py 146 i 50 96`, i = 0..9, primes 67108859, 67108837, 67108819, 67108777, 67108763, 67108757, 67108753, 67108747, 67108739, 67108729 (CRT in `sensed_146.json`).
* W_ref, W_rr: `run_wref.py 50 96 i`, `run_wrr.py 50 96 i`, i = 0..9, primes 16777213, 16777199, 16777183, 16777153, 16777141, 16777139, 16777127, 16777121, 16777099, 16777049.
* CRT range: `crt_bound.py` proves 2 R(v,f) < P for both tables, R = rooted c-nets (README §7).

## OEIS b-files in `data/`

Copies of the OEIS b-files the checks read, as downloaded in October 2026 (OEIS content, CC BY-SA 4.0).
`b000944.txt` is the OEIS file before the extension (a(1)-a(18)), not the new `b000944.txt` at the
top of this folder.

| file | source | terms | SHA-256 |
|---|---|---|---|
| `data/b000944.txt` | https://oeis.org/A000944/b000944.txt | 18 | `b67fea45a8c980ab4581d8c142132ecb7352d14dce4e6e3ef7bb6a897cab1482` |
| `data/b212438.txt` | https://oeis.org/A212438/b212438.txt | 196 | `a77fe823e03e6ae21c684a784fd94fdc08859f6b9fe1bb480157700aeab6b0d6` |
| `data/b239893.txt` | https://oeis.org/A239893/b239893.txt | 196 | `8745836dff2a156d0a9ee124ea6320679023b416837afabd89d4712c5c0d2a50` |
| `data/b000287.txt` | https://oeis.org/A000287/b000287.txt | 995 | `abe6016710df00bd50d018bb12e27a1dd152ae96bd74c0e96cc636fd69521392` |
| `data/b002840.txt` | https://oeis.org/A002840/b002840.txt | 30 | `6ed551d9ddf41198845214fd9b041e09fc1096dc8ff691a3503c7bff1b9a2f21` |
| `data/b006445.txt` | https://oeis.org/A006445/b006445.txt | 33 | `ff56eceefb5999aa03a8b4205595ce7aee7b0a2c2d432fdbbe8f76e6dda4ac96` |

## Residues of the a(50) run (`formula/`)

`combine_full.py 50 96 10` reads these and `sensed_146.json`; it rebuilds `Ufull_50_96.json` (published as `results/U_50_96.json`).

| file | SHA-256 |
|---|---|
| `formula/wref_50_96_16777049.npy` | `40e6f7c39f5fa89e1620046b65920dfce77abf53e5fd12f6116033adc4fdcde7` |
| `formula/wref_50_96_16777099.npy` | `bb951c6e606d4326514270e0dfb1e296c27e62d3c6c3871db8658851343f6604` |
| `formula/wref_50_96_16777121.npy` | `407f2956ee1c868c87420c80d726b6f43051efb607c23abdbeb0a824db514764` |
| `formula/wref_50_96_16777127.npy` | `d66ebdb11bd7439d566ea887641568e11bfbe549a05c23cde1eb9562c6a11918` |
| `formula/wref_50_96_16777139.npy` | `c9cf64cd63a38953f9868c72d5da8ee71789db4a743a0a81c00141eb5800cc98` |
| `formula/wref_50_96_16777141.npy` | `1b697612e143ec9ea43f4f432fd06bc331be200cb8ac5a335d2d07e87ac9f85d` |
| `formula/wref_50_96_16777153.npy` | `3f31887c9ee286a912e52fefdebd410f31a7e41db782587e8746a1a1e96e9c42` |
| `formula/wref_50_96_16777183.npy` | `cd3134e9e2b27f83fd4f01df1118093da0a8bc6e49732303707e11bc552c7c0b` |
| `formula/wref_50_96_16777199.npy` | `91e09296a92831e8bf5a8081095c1552a3aecdb4ec0daa606423ea36513e670f` |
| `formula/wref_50_96_16777213.npy` | `83eca31f45f4c462786a68d565875afd0354580169d44972d90f4b371c242c3d` |
| `formula/wrr_50_96_16777049.npy` | `d6dea15d5575ddea3ee1dc36ea770f0dce1e664a01320e0392316235c5159dff` |
| `formula/wrr_50_96_16777099.npy` | `6d8d7c16e974f2cfbd7eaae6b490f4c5bdb7947ff2fa424718a1506d2e210614` |
| `formula/wrr_50_96_16777121.npy` | `4a9b6b23dd865c6e34bc63d5d48318b7aeaab357d6ed2e551ba9ae562d7eb812` |
| `formula/wrr_50_96_16777127.npy` | `b6c5e0079548c91e3f64775c097557734a5aa66dfb603abe854328428da92399` |
| `formula/wrr_50_96_16777139.npy` | `8885f0dde697765461563f49d48484aded056c80104f9ddb6ec34fc22213841a` |
| `formula/wrr_50_96_16777141.npy` | `f1d0668cccccadae6758f9c34d7fe4ae7c2867815950aa5ca9d5e3d96569a5a1` |
| `formula/wrr_50_96_16777153.npy` | `40aadb1d047a3c60e0b30273fb66c41367e1e18d5ff78b974e6534376c256a5b` |
| `formula/wrr_50_96_16777183.npy` | `e9bf0ba5679c1bf3f98d8f4db0ec58f90f876638d1664aa89f13114dc01c1da9` |
| `formula/wrr_50_96_16777199.npy` | `c1196b025dfd92bafd946f484a4ad0cb97f8d22e7cfa3a6f0f0eeefb1ae586cb` |
| `formula/wrr_50_96_16777213.npy` | `acf252164f52cf4b859883b5f21e8e0fc025a4ad4cd68ff12090dd0bc8aa93d9` |
| `formula/sensed_146_67108729.npy` | `ef138af1d56167d554e97988530ae13bf982e431449db3ac7932252cb43a56e0` |
| `formula/sensed_146_67108739.npy` | `3a2810e0bf14450ada8900faf8048c8a065369e19b2e6ad506376a9159db4d1a` |
| `formula/sensed_146_67108747.npy` | `f1d8712a6f48162e0f233a4047cdfcff6317be124ef4f59bb2a01f3be4353d46` |
| `formula/sensed_146_67108753.npy` | `2253edfd0f056cac20dc404baf26640574770444f15e8daae2a1d05f2ceae159` |
| `formula/sensed_146_67108757.npy` | `f6d771bea674726bacd94cee908549d4f1c046a2653c554ab91c16301cd7ad61` |
| `formula/sensed_146_67108763.npy` | `88e9bbb122239328a6141f901483b563329ed3aacc832938a9fe0ef7f9e7a900` |
| `formula/sensed_146_67108777.npy` | `4978a59312e9e6469eae0b1ad2e37a54d9c577caf1971d2a1e9d840a9e9b8d97` |
| `formula/sensed_146_67108819.npy` | `ea75c0301c7808267d890c02b21ea969151f129297614c7f737c8bf5d5647664` |
| `formula/sensed_146_67108837.npy` | `dfc6db341200606e14d28dd54b308781f98318ee64ee3d9d9679e6bbc86a6ebd` |
| `formula/sensed_146_67108859.npy` | `6c11f2753488d1923a070c65b07522b1f8c5e910006692e9c13f6930365c4778` |

## Outputs

| file | SHA-256 |
|---|---|
| `b000944.txt` | `c3dddf12883d10f56d2a78298fe3ba65545024ab6525c44a8da9412de4d7bdee` |
| `results/U_50_96.json` | `a46528fc0513ba1a76fcb526212831adc08c9c6473817f43ceacc6b00347f25d` |
| `results/sensed_146.json` | `fed2ca8ce4bfeae30b2a346f802542b1d047ab69c21fbd1ef8e76cd2e7a6e7ae` |
| `results/a109_80.json` | `620a2284581b4a6933a07d5b4dcf5bab9d6fff5dbc65b46a7c9ca2fcd35ca4f6` |
