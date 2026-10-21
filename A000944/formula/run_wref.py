"""Compute W_ref mod one prime and cache it: python run_wref.py VX FY prime_index"""
import sys

from combine import wref_cached
from ms import no_throttle
from t3 import PRIMES24

if __name__ == '__main__':
    no_throttle()
    VX, FY, i = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    wref_cached(VX, FY, PRIMES24[i])
    print('done', VX, FY, PRIMES24[i])
