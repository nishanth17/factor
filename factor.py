import time
import math
import constants
import utils, primeSieve
import pollardRho, pollardPm1, ecm

small_primes = primeSieve.prime_sieve(constants.PRIME_THRESHOLD_BF)

def merge_factorizations(f1, f2):
	"""
	Merges prime factorizations of two numbers which are sorted in increasing order of 
	their prime factors into a larger one containing the prime factorization of their
	product -- similar to the merge step in mergesort. 
	"""
	if f1 == -1 or f2 == -1:
		# Factorization failed in this case
		return -1
	f = []
	i = j = 0
	while i < len(f1) and j < len(f2):
		if f1[i][0] < f2[j][0]:
			f.append(f1[i])
			i += 1
		elif f1[i][0] > f2[j][0]:
			f.append(f2[j])
			j += 1
		else:
			f.append((f1[i][0], f1[i][1] + f2[j][1]))
			i += 1
			j += 1
	if i < len(f1):
		f.extend(f1[i:])
	elif j < len(f2):
		f.extend(f2[j:])
	return f


def factorize_bf(n):
	"""
	Brute-forces small primes up to some pre-specified limit. 
	"""
	sn = int(math.sqrt(n))
	f = []
	for p in small_primes:
		if p > sn:
			if n > 1:
				f.append((n, 1))
				n = 1
			break
		i = 0
		while n % p == 0:
			n //= p
			i += 1
		if i > 0:
			f.append((p, i))
 			sn = int(math.sqrt(n))
	
	return f, n


def print_factoring_routine(n, routine_name):
	"""
	Prints factoring routine currently being used along with the number to be factored.  
	"""
	print "Factoring", str(n), "with", routine_name + "..."


# TODO: Incorporate Pollard (p-1) into this - ignoring it for now
def factorize(n, verbose = False, level = 3):
	"""
	Factorizes a specified integer or returns -1 if no factors can be found.
	"""
	if verbose: 
		if n != 1: 
			print "Factoring", str(n) + "..."
			print "Number of digits:", len(str(n))
	if n == 1:
		return []
	if utils.is_prime(n):
		if verbose:
			print str(n), "is prime!"
		return [(n, 1)]
	else:
		f, f1 = [], []
		if level > 2:
			# Try brute force for small prime factors
			if verbose: 
				print "Finding small prime factors..."
			f, n = factorize_bf(n)
			if verbose:
				if not f:
					print "Found no small prime factors... :("
				else:
					print "Prime factors found:", reduce(lambda x, y: x + y, [str(i[0]) + ", " for i in f])[:-2]

		
		if level > 1 and n <= constants.SIZE_THRESHOLD_RHO and n > 1:
			# Try Pollard rho
			if verbose:
				print_factoring_routine(n, constants.NAME_RHO)
			
			g = pollardRho.factorize_rho(n, verbose = verbose)
			if g != -1:
				if verbose:
					print "Found factor", str(g)
					f1 = merge_factorizations(factorize(g, verbose = verbose, level = 2), \
									factorize(n/g, verbose = verbose, level = 2))
					if f1 != -1:
						f.extend(f1)
		
		if level > 0 and (f1 == -1 or n > constants.SIZE_THRESHOLD_RHO) and n > 1:
			# If Pollard rho fails try ECM
			if verbose:
				print_factoring_routine(n, constants.NAME_ECM)

			g = ecm.factorize_ecm(n, verbose = verbose)
			if g != -1:
				if verbose:
					print "Found factor", str(g)
					f1 = merge_factorizations(factorize(g, verbose = verbose, level = 2), \
									factorize(n/g, verbose = verbose, level = 2))
					if f1 != -1:
						f.extend(f1)
					else:
						f = -1
		return f


def print_factorization(n, f):
	"""
	Prints a number as a product of the respective primes (and their exponents) in its prime 
	factorization.

	EXAMPLE:
	56 = 2^3 * 7^1
	"""
	if n == 1:
		return 1

	s = str(n) + " = "
	for i in xrange(len(f)-1):
		pf, exp = f[i][0], f[i][1]
		s += str(pf) + "^" + str(exp) + " * "
	
	s += str(f[-1][0]) + "^" + str(f[-1][1])
	return s


if __name__ == "__main__":
	while True:
		n = int(input("Enter number: "))
		print ""
		t = time.time()
		f = factorize(n, verbose = True)
		t1 = time.time()
		if f == -1:
			print "\n", n, "couldn't be factored :(\n"
		else:
			print "\n", print_factorization(n, f)
			print "\nTime:", t1 - t, "s\n"


		