"""Proposed DATA for each family entry, from results.json. Prints new DATA, new terms and length."""
import json
R = {k: {int(n): v for n, v in d.items()} for k, d in json.load(open('results.json')).items()}
UA, UL, US = R['unlabeled_all'], R['unlabeled_loopless'], R['unlabeled_simple']
LA, LL, LS = R['labeled_all'], R['labeled_loopless'], R['labeled_simple']


def tri(T, n0, k0, nmax=10):
    return [x for n in range(n0, nmax + 1) for x in T[n][k0:]]


# entry: (full sequence through n = 10, number of terms currently in DATA)
E = {
    'A002773': ([sum(US[n]) for n in range(11)], 10),
    'A058718': ([sum(UL[n]) for n in range(11)], 10),
    'A058730': (tri(US, 2, 2), 39),
    'A058716': (tri(UL, 0, 0), 45),
    'A058717': (tri(UL, 1, 1), 36),
    'A053534': (tri(UA, 0, 0), 55),
    'A058673': ([sum(LA[n]) for n in range(11)], 9),
    'A058669': (tri(LA, 0, 0), 47),
    'A058687': ([LA[n][3] for n in range(3, 11)], 6),
    'A058712': ([sum(LL[n]) for n in range(11)], 9),
    'A058710': (tri(LL, 0, 0), 45),
    'A058711': (tri(LL, 1, 1), 36),
    'A058721': ([sum(LS[n]) for n in range(2, 11)], 7),
    'A058720': (tri(LS, 2, 2), 28),
    'A058722': ([LS[n][4] for n in range(4, 11)], 5),
}
out = {}
for a, (seq, have) in E.items():
    data = ','.join(map(str, seq))
    out[a] = {'data': data, 'new': seq[have:], 'len': len(data)}
    print('%s len %3d new %s' % (a, len(data), seq[have:]))
json.dump(out, open('drafts.json', 'w'), indent=1)

# ---- b-files for the triangles whose rows 9-10 do not fit in DATA, and the submission draft ----
OFFSET = {'A058669': 0, 'A058710': 0, 'A058711': 1}
ROW = {'A058669': (LA, 0, 0), 'A058710': (LL, 0, 0), 'A058711': (LL, 1, 1)}
for a, off in OFFSET.items():
    seq = E[a][0]
    with open('b%s.txt' % a[1:], 'w', newline='\n') as f:
        for i, x in enumerate(seq):
            f.write('%d %d\n' % (i + off, x))
    T, n0, k0 = ROW[a]
    nine = tri(T, n0, k0, 9)
    out[a]['data'] = ','.join(map(str, nine))
    out[a]['new_in_data'] = nine[E[a][1]:]
    out[a]['len'] = len(out[a]['data'])
    out[a]['bfile'] = 'b%s.txt' % a[1:]
    out[a]['bfile_last'] = len(seq) - 1 + off
    assert out[a]['len'] <= 260
json.dump(out, open('drafts.json', 'w'), indent=1)
