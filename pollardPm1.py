# coding=utf-8

import math
import utils
import primeSieve

"""
This module contains an implementation of the two-stage variant of Pollard's 
p-1 algorithm.

This was adapted from a version of the same at StackExchange. 

References:
https://stackoverflow.com/questions/16424369/python-pollard-p-1-factorization

TODO:
Include explanation of algorithm. 
"""

MAX_B1 = 10**7
MAX_B2 = 10**9
MAX_D = 500

def compute_bounds(n):
	"""
	Computes Stage 1 and Stage 2 bounds for both Pollard p-1.
	"""
	log_q = math.log(pow(10, (len(str(n)) - 2) >> 1))
	t = int(math.ceil(math.exp(math.sqrt(0.5 * log_q * \
								math.log(log_q))) / 10) * 10)
	B1 = min(t, MAX_B1)
	B2 = min(B1 * 100, MAX_B2)
	return B1, B2


def factorize_pm1(n, verbose = False):
	B1, B2 = compute_bounds(n)
	if verbose:
		print "Bounds:", B1, B2
	
	primes_below_b1 = primeSieve.prime_sieve(B1)

	# ----- Stage 1 -----

	# Compute a large number which is B1-power-smooth. As in this implementation,
	# a usual choice for this number is the LCM of the integers below B1. 
	c = 2
	for p in primes_below_b1:
		pp = p
		while pp <= B1:
			c = pow(c, p, n)
			pp *= p

	g = utils.gcd(c-1, n)
	# If stage 1 is successful, return the non-trivial factor found. Else, go on
	# to stage 2. 
	if g != 1 and g != n:
		return g

	# ----- Stage 2 -----
	# NOTE: This stage only works if 'n' has exactly one prime factor between B1 and 
	# B2 (hence the name 'large-prime variant'). 
	d_cache = [-1] * (MAX_D+1)
	primes = primeSieve.segmented_sieve(B1+1, B2)
	p, temp_c = primes[0], c
	c, count = pow(c, p, n), 0

	for pos in xrange(1, len(primes)):
		q = primes[pos]
		# Use differences between successive primes and cache them
		d = q - p
		if d <= MAX_D:
			if d_cache[d] == -1:
				x = pow(temp_c, d, n)
				d_cache[d] = x
			else:
				x = d_cache[d]
		else:
			x = pow(temp_c, d, n)

		# Use modular multiplication instead of exponentiation to speed things up
		c, p = (c * x) % n, q
		count += 1

		# Accumulate products and compute GCD's periodically 
		if (count & 127) == 0:
			g = utils.gcd(c-1, n)
			# Return non-trivial factor if successful
			if g != 1 and g != n:
				return g

	g = utils.gcd(c-1, n)
	if g != 1 and g != n:
		return g
	else:
		return -1


if __name__ == "__main__":
	n = 3837523
	x = factorize_pm1(n, True)
	print n, x, n/x
