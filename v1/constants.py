# coding=utf-8

"""
File that contains constants tuned specifically for the factoring algorithms 
and the prime sieves. Tweakable if required. 
"""

# Prime sieve constants
SMALL_THRESHOLD = 60
ERAT_THRESHOLD = 35 * 10**5
ATKIN_THERSHOLD = 10**10
LOWER_SEG_SIZE = 65536
UPPER_SEG_SIZE = 2097152

# Pollard rho constants
PRIME_THRESHOLD_RHO = 500
SIZE_THRESHOLD_RHO = 10**20

# Pollard (p-1) constants
MAX_B1_PM1 = 10**8
MAX_B2_PM1 = 10**10
MAX_D_PM1 = 500

# ECM constants
MAX_CURVES_ECM = 10000
MAX_RND_ECM = 2**63
MAX_B1_ECM = 43 * 10**7
MAX_B2_ECM = 2 * 10**10

# General factorization constants
PRIME_THRESHOLD_BF = 25000

# Names of factoring routines for displaying purposes
NAME_ECM = "ECM"
NAME_RHO = "Pollard Rho"
NAME_PM1 = "Pollard p-1"