"""Write the 19 traced, symmetrized panels as text (panels.txt) and as JSON for the web page.
In panels.txt '#' is a printed (darned) mesh and '.' an open mesh, one panel per block."""
import glob, json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
src = sorted(glob.glob(os.path.join(HERE, '..', 'met1603_orig', 'panels', '*.npy')))
meta = json.load(open(os.path.join(HERE, '..', 'met1603_orig', 'panels.json')))
out, js = [], []
for f in src:
    name = os.path.basename(f)[:-4]; A = np.load(f).astype(int)
    page = int(name[1:3])
    out.append(f'# page {page} recto, panel {name[5:]}, {A.shape[0]} x {A.shape[0]}')
    out += [''.join('#' if v else '.' for v in row) for row in A]
    out.append('')
    js.append(dict(name=name, page=page, n=A.shape[0], rows=[''.join(str(v) for v in row) for row in A]))
open(os.path.join(HERE, '..', 'panels.txt'), 'w', newline='\n').write(
    '# Vinciolo, Les Secondes Oeuvres (Paris, 1603): the 19 fully symmetric square panels, traced from the\n'
    '# Metropolitan Museum of Art scans (39.30) and symmetrized by majority vote. # = printed mesh, . = open mesh.\n\n'
    + '\n'.join(out))
json.dump(js, open(os.path.join(HERE, '..', 'panels.json'), 'w'))
print(len(js), 'panels')
