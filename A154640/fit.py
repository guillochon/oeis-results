"""Fit rational g.f.s (shortest linear recurrence, verified on all remaining terms)."""
import sys, sympy as sp
x = sp.symbols('x')


def guess(a, maxd=None):
    maxd = maxd or (len(a) - 6) // 2
    for d in range(1, maxd + 1):
        for st in range(0, 6):
            seq = a[st:]
            if len(seq) < 2 * d + 4:
                continue
            A = sp.Matrix([[seq[i + j] for j in range(d)] for i in range(d)])
            if A.det() == 0:
                continue
            c = A.solve(sp.Matrix([seq[i + d] for i in range(d)]))
            if all(sum(c[j] * seq[i + j] for j in range(d)) == seq[i + d] for i in range(len(seq) - d)):
                den = 1 - sum(c[d - 1 - j] * x ** (j + 1) for j in range(d))
                num = sp.expand(sum(a[i] * x ** i for i in range(len(a))) * den)
                num = sum(num.coeff(x, i) * x ** i for i in range(d + st + 1))
                f = sp.factor(num / den)
                ser = sp.series(f, x, 0, len(a)).removeO()
                assert all(ser.coeff(x, i) == a[i] for i in range(len(a)))
                return d, st, len(a) - 2 * d - st, f
    return None


if __name__ == "__main__":
    for line in open(sys.argv[1]):
        k, v = line.split(' ', 1)
        print(k, guess(eval(v)))
