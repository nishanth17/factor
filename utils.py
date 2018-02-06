# coding=utf-8

import math
import random
import fractions

PRIME_THRESHOLD = 100000
MR_THRESHOLD = 10**36

def binary_search(x, arr, include_equal = False):
	"""
	Returns the index of the smallest element in an array which is larger
	than a specified element. This assumes that the array is sorted in 
	non-decreasing order. If the element is larger than the largest element 
	in the array, then the length of the array is returned instead. 

	Arguments:
		x (:int) - the element to be searched for
		arr (:int list) - the array sorted in non-decreasing order

	Returns:
		the position of the largest element in 'arr' greater than 'x'

	Examples:
		>>> binary_search(2, [0, 2, 3])
		>>> 2

		>>> binary_search(-1, [0, 2, 3])
		>>> 0

		>>> binary_search(99, [0, 2, 3])
		>>> 3
"""
	if x > arr[-1]:
		return len(arr)
	elif x < arr[0]:
		return 0

	l, r = 0, len(arr) - 1
	while l <= r:
		m = (l + r) >> 1
		if arr[m] == x:
			return m + 1 if not include_equal else m
		elif arr[m] < x:
			l = m + 1
		else:
			r = m - 1

	return l


def gcd(a, b):
	"""
	Returns the greatest common divisor (GCD) of two specified integers.

	Arguments:
		a (:int) - the first integer
		b (:int) - the second integer

	Reutrns:
		the GCD of 'a' and 'b'

	Examples:
		>>> gcd(1, 3)
		>>> 1

		>>> gcd(2, 4)
		>>> 2

		>>> gcd(10**8, 350)
		>>> 10
	"""
	return fractions.gcd(a, b)

def xgcd(a, b):
	"""
	Performs the Extended Euclidean algorithm to return the result of Bézout's 
	identity. 

	Arguments:
		a (:int) - the first integer
		b (:int) - the second integer

	Returns:
		'r' such that ar + bs = d where d = gcd(a, b)
	"""
	r, s = 0, 1
	while b != 0:
		c, d = divmod(a, b)
		r, s = s, r - c*s
		a, b = b, d
	return r


def is_prime_bf(n):
	"""
	Tests whether an integer is prime through brute force. A wheel (mod 6)
	is used to test potential candidates.

	Arguments:
		n (:int) - the integer to be tested

	Returns:
		True if 'n' is prime and False otherwise

	Examples:
		>>> is_prime_bf(20)
		>>> False

		>>> is_prime_bf(7)
		>>> True

		>>> is_prime_bf(9999)
		>>> False
	"""
	if n < 2: return False
	if n == 2 or n == 3: return True
	if not n & 1: return False
	if not n % 3: return False
	if n < 9: return True
	sqrt_n = int(math.sqrt(n)) + 1
	for i in range(5, sqrt_n, 6):
		if not n % i or not n % (i + 2): return False
	return True


def is_prime_fast(n, use_probabilistic = False, tolerance = 30):
	"""
	Tests whether a number is prime using a deterministic version of the Miller-
	Rabin primality test. Optionally tests whether the specified number is a 
	prime probabistically up to a given tolerance using the regular version of 
	the Miller-Rabin test. If the number is greater than 10^36, then all witnesses
	in the range [2, 2*log(n)*log(log(n))] are tested. However, this is conjectural
	and only heuristic evidence exists for it. To certify that a number is actually
	prime, one needs to test all witnesses in the range [2, 2*log(n)^2]. However, 
	this is generally quite slow. 

	Arguments:
		n (:int) - the integer to be tested
		use_probabilistic (:bool) - flag to indicate whether to use the regular 
		                		   version of the Miller-Rabin primality test
		tolerance (:int) - number of trials to be used to test primality

	Returns:
		True if 'n' is prime (or probably prime) and False otherwise

	Todo:
		Check for improved SPRP bases.

	References:
		- Francky from the PE Forums
		- https://miller-rabin.appspot.com/
		- https://en.wikipedia.org/wiki/Miller–Rabin_primality_test
	"""
	firstPrime = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, \
                  53, 59, 61, 67, 71]
    
    # Determine bases for deterministic Miller-Rabin test
	if n >= MR_THRESHOLD: 
		logn = math.log(n)
		if not use_probabilistic: 
			w = xrange(2, 2 * int(logn*log(logn)/log(2))) 
		else: 
			w = xrange(tolerance)
	elif n >= 1543267864443420616877677640751301: w = firstPrime[:20]
	elif n >= 564132928021909221014087501701: w = firstPrime[:18]
	elif n >= 59276361075595573263446330101: w = firstPrime[:16]
	elif n >= 6003094289670105800312596501: w = firstPrime[:15]
	elif n >= 3317044064679887385961981: w = firstPrime[:14]
	elif n >= 318665857834031151167461: w = firstPrime[:13]
	elif n >= 3825123056546413051: w = firstPrime[:12]
	#[2, 3, 5, 7, 11, 13, 17, 19, 23]
	elif n >= 341550071728321: w = firstPrime[:9]
	#[2, 3, 5, 7, 11, 13, 17]
	elif n >= 3474749660383: w = firstPrime[:7]
	elif n >= 2152302898749: w = firstPrime[:6]
	#[2, 3, 5, 7, 11, 13]
	elif n >= 4759123141: w = firstPrime[:5]
	#[2, 3, 5, 7, 11]
	elif n >= 9006403: w = [2, 7, 61]
	elif n >= 489997:
		# Some Fermat stuff
		if n&1 and n%3 and n%5 and n%7 and n%11 and n%13 and n%17 and n%19 \
		and n%23 and n%29 and n%31 and n%37 and n%41 and n%43 and n%47 \
		and n%53 and n%59 and n%61 and n%67 and n%71 and n%73 and n%79 \
		and n%83 and n%89 and n%97 and n%101:
			hn, nm1 = n >> 1, n - 1
			p = pow(2, hn, n)
			if p == 1 or p == nm1:
				p = pow(3, hn, n)
				if p == 1 or p == nm1:
					p = pow(5, hn, n)
					return p == 1 or p == nm1
		return False
	elif n >= 42799:
		return n&1 and n%3 and n%5 and n%7 and n%11 and n%13 and n%17 \
		and n%19 and n%23 and n%29 and n%31 and n%37 and n%41 and n%43 \
		and pow(2, n-1, n) == 1 and pow(5, n-1, n) == 1
	elif n >= 841:
		return n&1 and n%3 and n%5 and n%7 and n%11 and n%13 and n%17 \
		and n%19 and n%23 and n%29 and n%31 and n%37 and n%41 and n%43 \
		and n%47 and n%53 and n%59 and n%61 and n%67 and n%71 and n%73 \
		and n%79 and n%83 and n%89 and n%97 and n%101 and n%103 \
		and pow(2, n-1, n) == 1
	elif n >= 25:
		return n&1 and n%3 and n%5 and n%7 \
		and n%11 and n%13 and n%17 and n%19 and n%23
	elif n >= 4:
		return n&1 and n%3
	else:
		return n > 1
    
	if not (n&1 and n%3 and n%5 and n%7 and n%11 and n%13 and n%17 \
		   and n%19 and n%23 and n%29 and n%31 and n%37 and n%41 and n%43 \
		   and n%47 and n%53 and n%59 and n%61 and n%67 and n%71 and n%73 \
		   and n%79 and n%83 and n%89): return False
    
	# Miller-Rabin
	s = 0
	d = n - 1
	while not d & 1:
		d >>= 1
		s += 1
	for k in w:
		# Pick a random witness if probabilistic
		if use_probabilistic: 
			p = random.randint(2, n-2)
		else:
			p = k
		x = pow(p, d, n)
		if x == 1: continue
		for _ in xrange(s):
			if x+1 == n: break
			x = x*x % n
		else: return False
	return True


def is_prime(n, use_probabilistic = False, tolerance = 30):
	"""
	Tests whether a number is prime. The choice of test used depeneds on the size of 
	the specified number. Optionally tests whether the specified number is probably 
	prime up to a given tolerance using the regular version of the Miller-Rabin test. 

	Arguments:
		n (:int) - the integer to be tested
		use_probabilistic (:bool) - flag to indicate whether to use the regular 
		                		   version of the Miller-Rabin primality test
		tolerance (:int) - number of trials to be used to test primality

	Returns:
		True if 'n' is prime (or probably prime) and False otherwise

	Examples:
		>>> is_prime(20)
		>>> False

		>>> is_prime(7)
		>>> True

		>>> is_prime(9999)
		>>> False
	"""
	if n < PRIME_THRESHOLD: 
		return is_prime_bf(n)
	else: 
		if use_probabilistic:
			return is_prime_fast(n, use_probabilistic, tolerance)
		else:
			if n < MR_THRESHOLD:
				return is_prime_fast(n)
			else:
				return is_prime_fast(n, True, 40)

