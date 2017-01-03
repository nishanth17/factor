import math
import utils
import primeSieve
import pollardRho

SMALL_PRIME_THRESHOLD = 25000
SMALL_PRIMES = primeSieve.prime_sieve(SMALL_PRIME_THRESHOLD)

def merge_factorizations(f1, f2):
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


def factorize(n):
	if utils.is_prime(n):
		return [(n, 1)]
	else:
		# TODO: Try a bunch of factoring methods.
		pass