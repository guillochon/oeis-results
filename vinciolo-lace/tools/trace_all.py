import glob, sys, os
from multiprocessing import Pool
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
from trace import trace
def one(f):
    G, info, _ = trace(f)
    np.save(f[:-4] + '.npy', G)
    return f, info
if __name__ == '__main__':
    with Pool(6) as p:
        for f, info in sorted(p.map(one, sorted(glob.glob('p*r.jpg')))):
            print(os.path.basename(f), info)
