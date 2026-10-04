import utils
import random
import primeSieve
import constants

"""
This module contains an implementation of Brent's improvement to Pollard's Rho
alogrithm. 

TODO: Include explanation of algorithm. 
"""

small_primes = primeSieve.prime_sieve(constants.PRIME_THRESHOLD_RHO)

def factorize_rho(n, verbose = False):
    if n == 1 or utils.is_prime(n):
        return n

    # If no factor is found, return -1
    for i in range(len(small_primes) - 1, -1, -1):
        r, c, y = 1, small_primes[i], random.randint(1, n-1)
        if verbose:
            print "Trying offset:", c

        m, g, q, ys = random.randint(1, n-1), 1, 1, y
        min_val, k = 0, 0
        while g == 1:
            x, k = y, 0
            for j in range(r):
                y = y*y + c
                if y > n: y %= n
            while k < r and g == 1:
                ys, min_val = y, min(m, r-k)
                for j in range(min_val):
                    y = y*y + c
                    if y > n : y %= n
                    q = q * abs(x - y)
                    if q > n: q %= n
                g = utils.gcd(q, n)
                k += m
            r <<= 1
        
        if g == n:
            # If no factor found, try again.
            while True:
               ys = ys*ys + c
               if ys > n: ys %= n
               g = utils.gcd(abs(x-ys), n)
               if g > 1: 
                break
        
        if g != n:
            return g
        else:
            return -1
