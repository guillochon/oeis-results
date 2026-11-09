import json, numpy as np
def load():
    """19 panels as int arrays INCLUDING the frame (1 = printed/darned, 0 = open)."""
    return {p['name']: np.array([[int(c) for c in r] for r in p['rows']]) for p in json.load(open(__import__('os').path.join(__import__('os').path.dirname(__file__), '..', 'panels.json')))}
if __name__ == '__main__':
    for k, A in load().items():
        I = A[1:-1, 1:-1]; print(k, A.shape, 'frame ok' if A[0].all() and A[:,0].all() else 'NO FRAME', round(I.mean(), 3))
