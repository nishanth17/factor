import math
import utils
import primeSieve
import pollardRho
import constants

# Please don't change this
small_primes = primeSieve.prime_sieve(constants.PRIME_THRESHOLD_BF)

# Merges factorizations
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

# Returns divisors
def divisors_from_prime_factorization(f, sort = False):
    d = 1; r = []
    p = [0] * len(f)
    while True:
        r.append(d)
        i = 0
        while i < len(f) and p[i] == f[i][1]:
            p[i] = 0
            d //= f[i][0] ** f[i][1]
            i += 1
        if i >= len(f): break
        p[i] += 1
        d *= f[i][0]
    if sort: 
        return list(sorted(r))
    else:
        return r


def factorize(n):
	if utils.is_prime(n):
		return [(n, 1)]
	else:
		# TODO: Try a bunch of factoring methods.
		pass