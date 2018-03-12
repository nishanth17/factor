# coding=utf-8

""" 
This module has a bunch of prime sieves.

-> WHEELED SIEVE OF ERATOSTHENES
A sieve of Eratosthenes with a wheel mod 6.

-> SIEVE OF ATKIN
A segmented version of the sieve of Atkin as described in [1].
NOTE: This would probably be a lot more efficient with NumPy arrays but PyPy doesn't  
support NumPy as of yet. 

-> SEGMENTED SIEVE OF ERATOSTHENES
A segmented sieve of Eratosthenes with a wheel mod 2. The wheel mod 6 version of this is 
annoying as hell to implement and might be included in the future. 


BENCHMARKS:
Tests performed on a Macbook Pro (mid-2012) w/ a 2.6 GHz Intel Core i7 3720QM
processor and 8 GB RAM.

    BENCHMARKS |    10^6    |   10^7    |    10^8    |    10^9
    -------------------------------------------------------------
   Eratosthenes|    0.02s   |   0.32s   |    3.81s   |   93.99s
          Atkin|    0.06s   |   0.13s   |    0.72s   |    5.4s

REFRENCES:
[1] A.O.L Atkin, D.J.Bernstein; Prime Sieves using Binary Quadratic Forms; Mathematics  
    of Computation, 73-246: 1023-30
"""

import math
import time
import utils
import constants

# Sieve bits
segs = [[] for _ in xrange(60)]

# Primes under 60
under60 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59]

# delta's in the solutions to the congruences in algorithms 4.1, 4.2, 4.3
# in the paper
dAll = [1, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 49, 53, 59]

# All (d, f, g) where 4f^2 + g^2 = d (mod 60), d ≤ 60, f ≤ 15, g ≤ 30
DFG1 = [[1, 0, 1], [1, 0, 11], [1, 0, 19], \
			[1, 0, 29], [1, 2, 15], [1, 3, 5], [1, 3, 25], [1, 5, 9], \
			[1, 5, 21], [1, 7, 15], [1, 8, 15], [1, 10, 9], \
			[1, 10, 21], [1, 12, 5], [1, 12, 25], [1, 13, 15], \
			[13, 1, 3], [13, 1, 27], [13, 4, 3], [13, 4, 27], \
			[13, 6, 7], [13, 6, 13], [13, 6, 17], [13, 6, 23], \
			[13, 9, 7], [13, 9, 13], [13, 9, 17], [13, 9, 23], \
			[13, 11, 3], [13, 11, 27], [13, 14, 3], [13, 14, 27], \
			[17, 2, 1], [17, 2, 11], [17, 2, 19], [17, 2, 29], \
			[17, 7, 1], [17, 7, 11], [17, 7, 19], [17, 7, 29], \
			[17, 8, 1], [17, 8, 11], [17, 8, 19], [17, 8, 29], \
			[17, 13, 1], [17, 13, 11], [17, 13, 19], [17, 13, 29], \
			[29, 1, 5], [29, 1, 25], [29, 4, 5], [29, 4, 25], \
			[29, 5, 7], [29, 5, 13], [29, 5, 17], [29, 5, 23], \
			[29, 10, 7], [29, 10, 13], [29, 10, 17], [29, 10, 23], \
			[29, 11, 5], [29, 11, 25], [29, 14, 5], [29, 14, 25], \
			[37, 2, 9], [37, 2, 21], [37, 3, 1], [37, 3, 11], \
			[37, 3, 19], [37, 3, 29], [37, 7, 9], [37, 7, 21], \
			[37, 8, 9], [37, 8, 21], [37, 12, 1], [37, 12, 11], \
			[37, 12, 19], [37, 12, 29], [37, 13, 9], [37, 13, 21], \
			[41, 2, 5], [41, 2, 25], [41, 5, 1], [41, 5, 11], \
			[41, 5, 19], [41, 5, 29], [41, 7, 5], [41, 7, 25], \
			[41, 8, 5], [41, 8, 25], [41, 10, 1], [41, 10, 11], \
			[41, 10, 19], [41, 10, 29], [41, 13, 5], [41, 13, 25], \
			[49, 0, 7], [49, 0, 13], [49, 0, 17], [49, 0, 23], \
			[49, 1, 15], [49, 4, 15], [49, 5, 3], [49, 5, 27], \
			[49, 6, 5], [49, 6, 25], [49, 9, 5], [49, 9, 25], \
			[49, 10, 3], [49, 10, 27], [49, 11, 15], [49, 14, 15], \
			[53, 1, 7], [53, 1, 13], [53, 1, 17], [53, 1, 23], \
			[53, 4, 7], [53, 4, 13], [53, 4, 17], [53, 4, 23], \
			[53, 11, 7], [53, 11, 13], [53, 11, 17], [53, 11, 23], \
			[53, 14, 7], [53, 14, 13], [53, 14, 17], [53, 14, 23]]


# All (d, f, g) where 3f^2 + g^2 = d (mod 60), d ≤ 60, f ≤ 10, g ≤ 30
DFG2 = [[7, 1, 2], [7, 1, 8], [7, 1, 22], \
			[7, 1, 28], [7, 3, 10], [7, 3, 20], [7, 7, 10], \
			[7, 7, 20], [7, 9, 2], [7, 9, 8], [7, 9, 22], [7, 9, 28], \
			[19, 1, 4], [19, 1, 14], [19, 1, 16], [19, 1, 26], \
			[19, 5, 2], [19, 5, 8], [19, 5, 22], [19, 5, 28], \
			[19, 9, 4], [19, 9, 14], [19, 9, 16], [19, 9, 26], \
			[31, 3, 2], [31, 3, 8], [31, 3, 22], [31, 3, 28], \
			[31, 5, 4], [31, 5, 14], [31, 5, 16], [31, 5, 26], \
			[31, 7, 2], [31, 7, 8], [31, 7, 22], [31, 7, 28], \
			[43, 1, 10], [43, 1, 20], [43, 3, 4], [43, 3, 14], \
			[43, 3, 16], [43, 3, 26], [43, 7, 4], [43, 7, 14], \
			[43, 7, 16], [43, 7, 26], [43, 9, 10], [43, 9, 20]]


# All (d, f, g) where 3f^2 - g^2 = d (mod 60), d < 60, f ≤ 10, g ≤ 30
DFG3 = [[11, 0, 7], [11, 0, 13], [11, 0, 17], \
			[11, 0, 23], [11, 2, 1], [11, 2, 11], [11, 2, 19], \
			[11, 2, 29], [11, 3, 4], [11, 3, 14], [11, 3, 16], \
			[11, 3, 26], [11, 5, 2], [11, 5, 8], [11, 5, 22], \
			[11, 5, 28], [11, 7, 4], [11, 7, 14], [11, 7, 16], \
			[11, 7, 26], [11, 8, 1], [11, 8, 11], [11, 8, 19], \
			[11, 8, 29], [23, 1, 10], [23, 1, 20], [23, 2, 7], \
			[23, 2, 13], [23, 2, 17], [23, 2, 23], [23, 3, 2], \
			[23, 3, 8], [23, 3, 22], [23, 3, 28], [23, 4, 5], \
			[23, 4, 25], [23, 6, 5], [23, 6, 25], [23, 7, 2], \
			[23, 7, 8], [23, 7, 22], [23, 7, 28], [23, 8, 7], \
			[23, 8, 13], [23, 8, 17], [23, 8, 23], [23, 9, 10], \
			[23, 9, 20], [47, 1, 4], [47, 1, 14], [47, 1, 16], \
			[47, 1, 26], [47, 2, 5], [47, 2, 25], [47, 3, 10], \
			[47, 3, 20], [47, 4, 1], [47, 4, 11], [47, 4, 19], \
			[47, 4, 29], [47, 6, 1], [47, 6, 11], [47, 6, 19], \
			[47, 6, 29], [47, 7, 10], [47, 7, 20], [47, 8, 5], \
			[47, 8, 25], [47, 9, 4], [47, 9, 14], [47, 9, 16], \
			[47, 9, 26], [59, 0, 1], [59, 0, 11], [59, 0, 19], \
			[59, 0, 29], [59, 1, 2], [59, 1, 8], [59, 1, 22], \
			[59, 1, 28], [59, 4, 7], [59, 4, 13], [59, 4, 17], \
			[59, 4, 23], [59, 5, 4], [59, 5, 14], [59, 5, 16], \
			[59, 5, 26], [59, 6, 7], [59, 6, 13], [59, 6, 17], \
			[59, 6, 23], [59, 9, 2], [59, 9, 8], [59, 9, 22], \
			[59, 9, 28]]


def small_sieve(n):
	"""
	Returns the primes under a specified number with a modified sieve of Eratosthenes and a 
	wheel mod 6.

	Arguments:
		n (:int) - the number to list primes under

	Returns:
		the primes under 'n' in a list

	Examples:
		>>> small_sieve(9)
		>>> [2, 3, 5, 7]

		>>> small_sieve(30)
		>>> [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]

	References: 
		http://stackoverflow.com/questions/2068372/fastest-way-to-list-all-primes-below-n
		/3035188#3035188
	"""
	correction = (n % 6 > 1)
	n = {0: n, 1: n-1, 2: n+4, 3: n+3, 4: n+2, 5: n+1}[n % 6]
	sieve = [True] * (n/3)
	sieve[0] = False
	limit = int(math.sqrt(n))/3 + 1
	# Use a wheel (mod 6)
	for i in range(limit):
		if sieve[i]:
 			k = 3*i + 1 | 1
			sieve[((k*k)/3) :: (k << 1)] = \
					[False] * ((n/6 - (k*k)/6 - 1)/k + 1)
			sieve[(k * k + (k << 2) - \
					(k << 1) * (i & 1)) / 3 :: (k << 1)] = \
					[False] * ((n/6 - (k*k + (k << 2) - \
						2*k * (i & 1))/6 - 1)/k + 1)
	return [2, 3] + [3*i + 1 | 1 for i in xrange(1, n/3 - correction) if sieve[i]]


def enum1(d, f, g, L, B, segs):
	"""
	Alg 4.1: Given d ≤ 60, f ≤ 15, g ≤ 30 such that 4f^2 + g^2 = d (mod 60) find (x, y, k)
	with x > 0, y > 0, L ≤ k ≤ L + B, such that 4x^2 + y^2 = 60k + d and x = f + 15r, y = 
	g + 30s where r, s are integers.
	"""
	x, y0, temp = f, g, L+B
	k0 = (4*f*f + g*g - d) / 60
	while k0 < temp:
		k0 += x + x + 15
		x += 15

	while True:
		x -= 15
		k0 -= x + x + 15
		if x <= 0: 
			return
		while k0 < L:
			k0 += y0 + 15
			y0 += 30

		k, y = k0, y0
		while k < temp:
			segs[d][(k-L) >> 5] ^= 1 << ((k-L) & 31)
			k += y + 15
			y += 30


def enum2(d, f, g, L, B, segs):
	"""
	Alg 4.2: Given d ≤ 60, f ≤ 10, g ≤ 30 such that 3f^2 + g^2 = d (mod 60) find (x, y, k)
	with x > 0, y > 0, L ≤ k ≤ L + B, such that 3x^2 + y^2 = 60k + d and x = f + 10r, y = 
	g + 30s where r, s are integers.
	"""
	x, y0, temp = f, g, L+B
	k0 = (3*f*f + g*g - d) / 60
	while k0 < temp:
		k0 += x + 5
		x += 10

	while True:
		x -= 10
		k0 -= x + 5
		if x <= 0: 
			return
		while k0 < L:
			k0 += y0 + 15
			y0 += 30

		k, y = k0, y0
		while k < temp:
			segs[d][(k-L) >> 5] ^= 1 << ((k-L) & 31)

			k += y + 15
			y += 30


def enum3(d, f, g, L, B, segs):
	"""
	Alg 4.3: Given d < 60, f ≤ 10, g ≤ 30 such that 3f^2 - g^2 = d (mod 60) find (x, y, k)
	with x > 0, y > 0, L ≤ k ≤ L + B, such that 3x^2 - y^2 = 60k + d and x = f + 10r, y = 
	g + 30s where r, s are integers.
	"""
	x, y0, temp = f, g, L+B
	k0 = (3*f*f - g*g - d) / 60

	while True:
		while k0 >= temp:
			if x <= y0: 
				return
			k0 -= y0 + 15
			y0 += 30

		k, y = k0, y0
		while k >= L and y < x:
			segs[d][(k-L) >> 5] ^= 1 << ((k-L) & 31)
			k -= y + 15
			y += 30
		
		k0 += x + 5
		x += 10


def sieve_of_atkin(n):
	"""
	Returns the primes under a specified number with a segmented sieve of Atkin.

	Arguments:
		n (:int) - the number to list primes under

	Returns:
		the primes under 'n' in a list
	"""
	sqrt_n, u, r = int(math.sqrt(n)), n + 32, 17
	B, lu = 60 * sqrt_n, math.log(u)
	primes = small_sieve(sqrt_n)
	ret = under60 + [0] * int(u/lu + u/(lu*lu) * 1.5 - r)
	for d in dAll:
		segs[d] = [0] * ((B >> 5) + 1)

	# Do computations in segments of size 60√n
	lim = n/60 + 1
	for L in xrange(1, lim, B):
		for d in dAll:
			for k in xrange(len(segs[d])):
				segs[d][k] = 0

		# Sieve off the primes (i.e. solutions to the various quadratic
		# Diophantine equations)
		lim2 = 60 * (L+B)
		for d,f,g in DFG1:
			enum1(d, f, g, L, B, segs)
		for d,f,g in DFG2:
			enum2(d, f, g, L, B, segs)
		for d,f,g in DFG3:
			enum3(d, f, g, L, B, segs)

		# Sieve off non-squarefree numbers
		for p in primes:
			p2 = p * p
			if p2 > lim2: 
				break
			if p >= 7:
				b = -utils.xgcd(p2, 60)
				if b < 0: b += p2
				for d in dAll:
					x = b * (60*L + d) % p2
					while x < B:
						segs[d][x >> 5] &= ~(1 << (x & 31))
						x += p2

		# Compute primes
		for j in xrange((B >> 5) + 1):
			for x in xrange(32):
				k = 60 * (L + x + (j << 5))
				for d in dAll:
					if k + d > n:
						return ret[:r]
					# If a_k = 1, 60k + d is a prime
					if ((segs[d][j] << 31 - x) & 0xFFFFFFFF) >= 0x80000000:
						ret[r] = 60*k + d
						r += 1

def prime_sieve(n):
	"""
	Returns the primes below a specified number with the choice of prime sieve depending on the 
	size of the number. 

	Arguments:
		n (:int) - the number to list primes under

	Returns:
		the primes under 'n' in a list

	Examples:
		>>> prime_sieve(9)
		>>> [2, 3, 5, 7]

		>>> len(prime_sieve(10**9))
		>>> 50847534
	"""
	if n <= constants.SMALL_THRESHOLD:
		return under60[:utils.binary_search(n, under60)]
	elif n <= constants.ERAT_THRESHOLD:
		return small_sieve(n)
	elif n <= constants.ATKIN_THERSHOLD:
		return sieve_of_atkin(n)
	else:
		return segmented_sieve(2, n)


def segmented_sieve(lo, hi):
	"""
	Returns the primes between two specified numbers using a segmented sieve of Eratosthenes. 
	Optionally, one may specify the size of the segment to be used. If not specified, the segment 
	size used defaults to the square root of the difference between the two specified numbers.

	NOTE: A small segment size results in low memory usage but results in a large computation time.
	There seems to be an optimal segment size but I can't really figure out what it is.

	Arguments:
		lo (:int) - the lower bound of the interval
		hi (:int) - the upper bound of the interval

	Returns:
		the primes in the interval [lo, hi] in a list
	"""
	if hi < lo: return []
	max_prime, pos = int(math.sqrt(hi)), 0
	base_primes = prime_sieve(max_prime)
	primes = [0] * int(math.ceil(1.5 * hi/math.log(hi)) - math.floor(1.5 * lo/math.log(lo)))

	# Include primes below √hi if necessary
	if lo < max_prime:
		lo_pos = utils.binary_search(lo, base_primes, include_equal = True)
		for k in xrange(lo_pos, len(base_primes)):
			primes[pos] = base_primes[k]
			pos += 1
		lo = max_prime

	# Compute segment size 
	delta = constants.UPPER_SEG_SIZE if hi - lo >= constants.UPPER_SEG_SIZE else constants.LOWER_SEG_SIZE

	l1, l = len(base_primes), (delta >> 4) + 1
	int_size, sieve = l << 3, bytearray([0x0] * l)
	lo_1, hi_1 = lo, lo + delta
	
	# Compute stuff in segments
	while lo_1 <= hi:
		# Re-zero sieve bits if necessary
		if lo_1 != lo:
			for i in range(l):
				sieve[i] = 0

		if (lo_1 & 1) == 0: 
			lo_1 += 1

		# Sieve off primes
		for i in xrange(1, l1):
			p = base_primes[i]
			k = (p - (lo_1 % p)) % p
			if (k & 1) == 1: 
				k += p
			k >>= 1
			while k < int_size:
				sieve[k >> 3] |= 1 << (k & 7)
				k += p

		# Compute primes and put them in the prime list
		end = min(hi_1, hi) + 1
		for n in range(lo_1, end, 2):
			d = n - lo_1
			if ((sieve[d >> 4] >> ((d >> 1) & 0x7)) & 0x1) == 0x0:
				primes[pos] = n
				pos += 1

		# Update segment boundaries
		lo_1 = hi_1 + 1
		hi_1 = lo_1 + delta

	return primes[:pos]

if __name__ == "__main__":
	N = 10**8
	t1 = time.time()
	x = segmented_sieve(2, N)
	t2 = time.time()
	print "Num primes:", len(x)
	print "Time:", (t2-t1), "s"

