# coding=utf-8

import math
import utils
import random
import primeSieve
import constants
from decimal import Decimal

"""
This module contains an implementation of a two-stage version Lenstra's elliptical 
curve factorization method (ECM) with the usual stage 1 and stage 2 optimizations. 
This implementation uses Suyama's paramterization to generate curves in Montgomery
form and is inversionless.
"""

RESOLUTION = 40

def compute_bounds(n):
	"""
	Computes Stage 1 and Stage 2 bounds for both ECM. This almost  coincides with GMP-ECM's 
	bounds for the same but are clipped because of prime sieve limitations. 
	Reference:
	http://www.mersennewiki.org/index.php/Elliptic_Curve_Method	
	"""
	log_n = len(str(n))
	if log_n <= 30: 
		B1, B2 = 2000, 147396
	elif log_n <= 40: 
		B1, B2 = 11000, 1873422
	elif log_n <= 50: 
		B1, B2 = 50000, 12746592
	elif log_n <= 60: 
		B1, B2 = 250000, 128992510
	elif log_n <= 70: 
		B1, B2 = 1000000, 1045563762
	elif log_n <= 80:
		B1, B2 = 3000000, 5706890290
	else: 
		# Anything greater and my computer runs out of memory -- prolly need to fix this
		B1, B2 = constants.MAX_B1_ECM, constants.MAX_B2_ECM
	return B1, B2


def point_add(px, pz, qx, qz, rx, rz, n):
	"""
	Adds two specified P and Q points (in Montgomery form) in E(Z\nZ). Assumes R = P - Q.
	"""
	u = (px-pz) * (qx+qz)
	v = (px+pz) * (qx-qz)
	upv, umv = u+v, u-v
	x = (rz * upv * upv)
	if x >= n:
		x %= n
	z = rx * umv * umv
	if z >= n:
		z %= n
	return x, z


def point_double(px, pz, n, a24):
	"""
	Doubles a point P (in Montgomery form) in E(Z\nZ).
	"""
	u, v = px+pz, px-pz
	u2, v2 = u*u, v*v
	t = u2 - v2
	x = (u2 * v2) 
	if x >= n:
		x %= n
	z = (t * (v2 + a24*t))
	if z >= n:
		z %= n
	return x, z


def scalar_multiply(k, px, pz, n, a24):
	"""
	Multiplies a specified point P (in Montgomery form) by a specified scalar in E(Z\nZ).
	"""
	sk = bin(k)
	lk = len(sk)
	qx, qz = px, pz
	rx, rz = point_double(px, pz, n, a24)

	for i in xrange(3, lk):
		if sk[i] == '1':
			qx, qz = point_add(rx, rz, qx, qz, px, pz, n)
			rx, rz = point_double(rx, rz, n, a24)
		else:
			rx, rz = point_add(qx, qz, rx, rz, px, pz, n)
			qx, qz = point_double(qx, qz, n, a24)	

	return qx, qz

###########################################################

ADD_COST = 6
DUP_COST = 5

def lucas_cost(k, v):
	d = k
	r = int(Decimal(d) * Decimal(v) + Decimal(0.5))
	if r >= k:
		return ADD_COST * k

	d, e, c = k - r, 2*r - k, DUP_COST + ADD_COST
	while d != e:
		# Want d >= e so swap if d < e
		if d < e:
			d, e = e, d

		# Condition 1
		if 4*d <= 5*e and (d + e) % 3 == 0:
			d, e = (2*d - e) / 3, (2*e - d) / 3
			c += 3 * ADD_COST
		# Condition 2
		elif 4*d <= 5*e and (d - e) % 6  == 0:
			d = (d - e) / 2
			c += ADD_COST + DUP_COST
		# Condition 3
		elif d <= 4*e:
			d -= e
			c += ADD_COST
		# Condition 4
		elif (d + e) % 2 == 0:
			d = (d - e) / 2
			c += ADD_COST + DUP_COST
		# Condition 5
		elif d % 2 == 0:
			d /= 2
			c += ADD_COST + DUP_COST
		# Condition 6
		elif d % 3 == 0:
			d = d/3 - e
			c += 3*ADD_COST + DUP_COST
		# Condition 7
		elif (d + e) % 3 == 0:
			d = (d - 2*e) / 3
			c += 3*ADD_COST + DUP_COST
		# Condition 8
		elif (d - e) % 3 == 0:
			d = (d - e) / 3
			c += 3*ADD_COST + DUP_COST
		# Condition 9
		else:
			e /= 2
			c += ADD_COST + DUP_COST

	return c


def multiply_prac(k, px, pz, n, a24):
	ax, bx, cx, tx, t2x = px, 0, 0, 0, 0
	az, bz, cz, tz, t2z = pz, 0, 0, 0, 0
	v = [0.61803398874989485, 0.5801787282954641, 0.6179144065288179 , 0.6180796684698958]

	# Find best value of v
	r, i = lucas_cost(k, v[0]), 0
	for d in xrange(len(v)):
		e = lucas_cost(k, v[d])
		if e < r:
			r, i = e, d

	r = int(Decimal(k) * Decimal(v[i]) + Decimal(0.5))
	d, e = k - r, 2*r - k
	bx, bz, cx, cz = ax, az, ax, az
	ax, az = point_double(ax, az, n, a24)
	
	while d != e:
		# Want d >= e so swap if d < e
		if d < e:
			d, e = e, d
			ax, az, bx, bz = bx, bz, ax, az

		# Condition 1
		if 4*d <= 5*e and (d + e) % 3 == 0:
			d, e = (2*d - e) / 3, (2*e - d) / 3
			tx, tz = point_add(ax, az, bx, bz, cx, cz, n)
			t2x, t2z = point_add(tx, tz, ax, az, bx, bz, n)
			bx, bz = point_add(bx, bz, tx, tz, ax, az, n)
			ax, az, t2x, t2z = t2x, t2z, ax, az
		# Condition 2
		elif 4*d <= 5*e and (d - e) % 6  == 0:
			d = (d - e) / 2
			bx, bz = point_add(ax, az, bx, bz, cx, cz, n)
			ax, az = point_double(ax, az, n, a24)
		# Condition 3 
		elif d <= 4*e:
			d -= e
			# tx, tz = point_add(bx, bz, ax, az, cx, cz, n)
			# bx, tx, cx = tx, cx, bx
			# bz, tz, cz = tz, cz, bz
			cx, cz = point_add(bx, bz, ax, az, cx, cz, n)
			bx, bz, cx, cz = cx, cz, bx, bz
		# Condition 4
		elif (d + e) % 2 == 0:
			d = (d - e) / 2
			bx, bz = point_add(bx, bz, ax, az, cx, cz, n)
			ax, az = point_double(ax, az, n, a24)
		# Condition 5
		elif d % 2 == 0:
			d /= 2
			cx, cz = point_add(cx, cz, ax, az, bx, bz, n)
			ax, az = point_double(ax, az, n, a24)
		# Condition 6
		elif d % 3 == 0:
			d = d/3 - e
			tx, tz = point_double(ax, az, n, a24)
			t2x, t2z = point_add(ax, az, bx, bz, cx, cz, n)
			ax, az = point_add(tx, tz, ax, az, ax, az, n)
			# tx, tz = point_add(tx, tz, t2x, t2z, cx, cz, n)
			# cx, bx, tx = bx, tx, cx
			# cz, bz, tz = bz, tz, cz
			cx, cz = point_add(tx, tz, t2x, t2z, cx, cz, n)
			bx, bz, cx, cz = cx, cz, bx, bz
		# Condition 7
		elif (d + e) % 3 == 0:
			d = (d - 2*e) / 3
			tx, tz = point_add(ax, az, bx, bz, cx, cz, n)
			bx, bz = point_add(tx, tz, ax, az, bx, bz, n)
			tx, tz = point_double(ax, az, n, a24)
			# TODO: Check order of a and t here
			ax, az = point_add(ax, az, tx, tz, ax, az, n)
		# Condition 8
		elif (d - e) % 3 == 0:
			d = (d - e) / 3
			tx, tz = point_add(ax, az, bx, bz, cx, cz, n)
			# TODO: Check whether c = f(a, c, b) or c = f(c, a, b)
			cx, cz = point_add(cx, cz, ax, az, bx, bz, n)
			bx, bz, tx, tz = tx, tz, bx, bz
			tx, tz = point_double(ax, az, n, a24)
			# TODO: Check order of a and t here
			ax, az = point_add(ax, az, tx, tz, ax, az, n)
		# Condition 9
		else:
			e /= 2
			cx, cz = point_add(cx, cz, bx, bz, ax, az, n)
			bx, bz = point_double(bx, bz, n, a24)
	
	x, z = point_add(ax, az, bx, bz, cx, cz, n)
	return x, z


###########################################################


def factorize_ecm(n, verbose = False):
	"""
	ECM algorithm
	"""
	if n == 1 or utils.is_prime(n):
		return n
        
	B1, B2 = compute_bounds(n)
	if verbose:
		print "Number of digits:", len(str(n))
		print "Bounds:", B1, B2

	D = int(math.sqrt(B2))
	beta = [0] * (D+1)
	S = [0] * (2*D + 2)

	# ----- Stage 1 and Stage 2 precomputations -----
	curves, log_B1 = 0, math.log(B1)

	if verbose: print "Sieving primes..."
	primes = primeSieve.prime_sieve(B2)

	num_primes = len(primes)
	idx_B1 = utils.binary_search(B1, primes)
	
	# Compute a B1-powersmooth integer 'k'
	k = 1
	for i in xrange(idx_B1):
		p = primes[i]
		k = k * pow(p, int(log_B1/math.log(p)))

	g = 1
	while (g == 1 or g == n) and curves <= constants.MAX_CURVES_ECM:
		curves += 1
		sigma = random.randint(6, constants.MAX_RND_ECM)
		if verbose and curves % RESOLUTION == 0: 
			print "Tried", curves, "random curves..."

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
			print "Stage 1 found factor!"
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

	# No non-trivial factor found, return -1
	if curves > constants.MAX_CURVES_ECM:
		return -1
	else:
		print "Stage 2 found factor!"
		return g