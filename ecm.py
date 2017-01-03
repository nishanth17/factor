import math
import utils
import random
import primeSieve

"""
This module contains an implementation of a two-stage version Lenstra's elliptical 
curve factorization method (ECM) with the usual stage 1 and stage 2 optimizations. 
"""

MAX_B1 = 2 * 10**7
MAX_B2 = 2 * 10**9
MAX_CURVES = 5000
MAX_RND = (1 << 32) - 1

def compute_bounds(n):
	"""
	Computes Stage 1 and Stage 2 bounds.
	"""
	log_q = math.log(pow(10, (len(str(n)) - 2) >> 1))
	t = int(math.ceil(math.exp(math.sqrt(0.5 * log_q * \
								math.log(log_q))) / 10) * 10)
	B1 = min(t, MAX_B1)
	B2 = min(B1 * 100, MAX_B2)
	return B1, B2


def point_add(px, pz, qx, qz, rx, rz, n):
	"""
	Adds three specified points (in Montgomery form) in E(Z\nZ).
	"""
	u = (px-pz) * (qx+qz)
	v = (px+pz) * (qx-qz)
	upv, umv = u+v, u-v
	x = (rz * upv * upv) % n
	z = (rx * umv * umv) % n
	return x, z


def point_double(px, pz, n, a24):
	"""
	Doubles a point (in Montgomery form) in E(Z\nZ).
	"""
	u, v = px+pz, px-pz
	u2, v2 = u*u, v*v
	t = u2 - v2
	x = (u2 * v2) % n
	z = (t * (v2 + a24*t)) % n
	return x, z


def scalar_multiply(k, px, pz, n, a24):
	"""
	Multiplies a specified point (in Montgomery form) by a specified scalar in E(Z\nZ).
	"""
	sk = bin(k)
	lk = len(sk)
	qx, qz = px, pz
	rx, rz = point_double(px, pz, n, a24)

	for i in xrange(3, lk):
		if sk[i] == 1:
			qx, qz = point_add(rx, rz, qx, qz, px, pz, n)
			rx, rz = point_double(rx, rz, n, a24)
		else:
			rx, rz = point_add(qx, qz, rx, rz, px, pz, n)
			qx, qz = point_double(qx, qz, n, a24)	

	return qx, qz


def factorize_ecm(n, verbose = False):
	"""
	ECM stuff.
	"""
	B1, B2 = compute_bounds(n)
	if verbose:
		print "Bounds:", B1, B2

	D = int(math.sqrt(B2))
	beta = [0] * (D+1)
	S = [0] * (2*D + 2)

	# ----- Stage 1 and Stage 2 precomputations -----
	curves, log_B1 = 0, math.log(B1)
	primes = primeSieve.prime_sieve(B2)
	num_primes = len(primes)
	idx_B1 = utils.binary_search(B1, primes)
	
	# Compute a B1-powersmooth integer 'k'
	k = 1
	for i in xrange(idx_B1):
		p = primes[i]
		k = k * pow(p, int(log_B1/math.log(p)))

	g = 1
	while (g == 1 or g == n) and curves <= MAX_CURVES:
		curves += 1
		sigma = random.randint(6, MAX_RND)
		if verbose:
			print "Curve", str(curves) + ":", sigma

		# Generate a new random curve in Montgomery form with Suyama's parametrization
		u = ((sigma * sigma) - 5) % n
		v = (4 * sigma) % n
		vmu = v - u
		A = ((vmu*vmu*vmu) * (3*u + v) / (4*u*u*u*v) - 2) % n
		a24 = (A+2) / 4

		# ----- Stage 1 -----
		px, pz = ((u*u*u) / (v*v*v)) % n, 1
		qx, qz = scalar_multiply(k, px, pz, n, a24)
		g = utils.gcd(n, qz)

		# If stage 1 is successful, return a non-trivial factor else
		# move on to stage 2
		if g != 1 and g != n:
			return g

		# ----- Stage 2 -----
		S[1], S[2] = point_double(qx, qz, n, a24)
		S[3], S[4] = point_double(S[1], S[2], n, a24)
		beta[1] = (S[1] * S[2]) % n
		beta[2] = (S[3] * S[4]) % n
		for d in xrange(3, D+1):
			d2 = 2 * d
			S[d2-1], S[d2] = point_add(S[d2-3], S[d2-2], S[1], S[2], S[d2-5], S[d2-4], n)
			beta[d] = (S[d2-1] * S[d2]) % n

		g, B = 1, B1 - 1
		rx, rz = scalar_multiply(B, qx, qz, n, a24)
		tx, tz = scalar_multiply(B - 2*D, qx, qz, n, a24)
		
		q, step = idx_B1, 2*D
		for r in xrange(B, B2, step):
			alpha, limit = (rx * rz) % n, r + step
			
			while q < num_primes and primes[q] <= limit:
				d = (primes[q] - r) / 2
				f = (rx - S[2*d-1]) * (rz + S[2*d]) - alpha + beta[d]
				g = (g * f) % n
				q += 1

			trx, trz = rx, rz
			rx, rz = point_add(rx, rz, S[2*D-1], S[2*D], tx, tz, n)
			tx, tz = trx, trz

		g = utils.gcd(n, g)

	# No non-trivial factor found, return 0
	if curves > MAX_CURVES:
		return 0
	else:
		return g