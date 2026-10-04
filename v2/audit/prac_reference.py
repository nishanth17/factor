"""Audit prototype: guarded PRAC chain selection using the original point updates.
Not a complete elliptic-curve group API; see report for exceptional-point limits.
"""

from math import gcd

ADD_COST = 6
DUP_COST = 5

def point_add(px, pz, qx, qz, rx, rz, n):
    """
        Adds two specified P and Q points (in Montgomery form) in E(Z
Z). Assumes R = P - Q.
        """
    u = (px - pz) * (qx + qz)
    v = (px + pz) * (qx - qz)
    (upv, umv) = (u + v, u - v)
    x = rz * upv * upv
    if x >= n:
        x %= n
    z = rx * umv * umv
    if z >= n:
        z %= n
    return (x, z)

def point_double(px, pz, n, a24):
    """
        Doubles a point P (in Montgomery form) in E(Z
Z).
        """
    (u, v) = (px + pz, px - pz)
    (u2, v2) = (u * u, v * v)
    t = u2 - v2
    x = u2 * v2
    if x >= n:
        x %= n
    z = t * (v2 + a24 * t)
    if z >= n:
        z %= n
    return (x, z)

def scalar_multiply(k, px, pz, n, a24):
    """
        Multiplies a specified point P (in Montgomery form) by a specified scalar in E(Z
Z).
        """
    if k < 0:
        raise ValueError('negative scalar')
    if k == 0:
        return (1, 0)
    sk = bin(k)
    lk = len(sk)
    (qx, qz) = (px, pz)
    (rx, rz) = point_double(px, pz, n, a24)
    for i in range(3, lk):
        if sk[i] == '1':
            (qx, qz) = point_add(rx, rz, qx, qz, px, pz, n)
            (rx, rz) = point_double(rx, rz, n, a24)
        else:
            (rx, rz) = point_add(qx, qz, rx, rz, px, pz, n)
            (qx, qz) = point_double(qx, qz, n, a24)
    return (qx, qz)

def _chain_cost(k, r):
    (d, e, c) = (k - r, 2 * r - k, DUP_COST + ADD_COST)
    while d != e:
        if d < e:
            (d, e) = (e, d)
        if 4 * d <= 5 * e and (d + e) % 3 == 0:
            (d, e) = ((2 * d - e) // 3, (2 * e - d) // 3)
            c += 3 * ADD_COST
        elif 4 * d <= 5 * e and (d - e) % 6 == 0:
            d = (d - e) // 2
            c += ADD_COST + DUP_COST
        elif d <= 4 * e:
            d -= e
            c += ADD_COST
        elif (d + e) % 2 == 0:
            d = (d - e) // 2
            c += ADD_COST + DUP_COST
        elif d % 2 == 0:
            d //= 2
            c += ADD_COST + DUP_COST
        elif d % 3 == 0:
            d = d // 3 - e
            c += 3 * ADD_COST + DUP_COST
        elif (d + e) % 3 == 0:
            d = (d - 2 * e) // 3
            c += 3 * ADD_COST + DUP_COST
        elif (d - e) % 3 == 0:
            d = (d - e) // 3
            c += 3 * ADD_COST + DUP_COST
        else:
            e //= 2
            c += ADD_COST + DUP_COST
    return c

RATIOS = (61803398874989485, 58017872829546410, 61791440652881790, 61807966846989580)
DENOM = 10**17

def select_chain(k):
    candidates = set()
    for numerator in RATIOS:
        r = (k*numerator + DENOM//2)//DENOM
        for candidate in (r-1, r, r+1):
            if 0 < 2*candidate-k and candidate < k and gcd(k,candidate) == 1:
                candidates.add(candidate)
    if not candidates:
        return None
    return min(candidates, key=lambda r: (_chain_cost(k,r),r))


def multiply_prac(k, px, pz, n, a24):
    if k < 0:
        raise ValueError('negative scalar')
    if k == 0:
        return (1, 0)
    if k == 1:
        return (px, pz)
    if k == 2:
        return point_double(px, pz, n, a24)
    r = select_chain(k)
    if r is None:
        return scalar_multiply(k, px, pz, n, a24)
    (ax, bx, cx, tx, t2x) = (px, 0, 0, 0, 0)
    (az, bz, cz, tz, t2z) = (pz, 0, 0, 0, 0)
    (d, e) = (k - r, 2 * r - k)
    (bx, bz, cx, cz) = (ax, az, ax, az)
    (ax, az) = point_double(ax, az, n, a24)
    while d != e:
        if d < e:
            (d, e) = (e, d)
            (ax, az, bx, bz) = (bx, bz, ax, az)
        if 4 * d <= 5 * e and (d + e) % 3 == 0:
            (d, e) = ((2 * d - e) // 3, (2 * e - d) // 3)
            (tx, tz) = point_add(ax, az, bx, bz, cx, cz, n)
            (t2x, t2z) = point_add(tx, tz, ax, az, bx, bz, n)
            (bx, bz) = point_add(bx, bz, tx, tz, ax, az, n)
            (ax, az, t2x, t2z) = (t2x, t2z, ax, az)
        elif 4 * d <= 5 * e and (d - e) % 6 == 0:
            d = (d - e) // 2
            (bx, bz) = point_add(ax, az, bx, bz, cx, cz, n)
            (ax, az) = point_double(ax, az, n, a24)
        elif d <= 4 * e:
            d -= e
            (cx, cz) = point_add(bx, bz, ax, az, cx, cz, n)
            (bx, bz, cx, cz) = (cx, cz, bx, bz)
        elif (d + e) % 2 == 0:
            d = (d - e) // 2
            (bx, bz) = point_add(bx, bz, ax, az, cx, cz, n)
            (ax, az) = point_double(ax, az, n, a24)
        elif d % 2 == 0:
            d //= 2
            (cx, cz) = point_add(cx, cz, ax, az, bx, bz, n)
            (ax, az) = point_double(ax, az, n, a24)
        elif d % 3 == 0:
            d = d // 3 - e
            (tx, tz) = point_double(ax, az, n, a24)
            (t2x, t2z) = point_add(ax, az, bx, bz, cx, cz, n)
            (ax, az) = point_add(tx, tz, ax, az, ax, az, n)
            (cx, cz) = point_add(tx, tz, t2x, t2z, cx, cz, n)
            (bx, bz, cx, cz) = (cx, cz, bx, bz)
        elif (d + e) % 3 == 0:
            d = (d - 2 * e) // 3
            (tx, tz) = point_add(ax, az, bx, bz, cx, cz, n)
            (bx, bz) = point_add(tx, tz, ax, az, bx, bz, n)
            (tx, tz) = point_double(ax, az, n, a24)
            (ax, az) = point_add(ax, az, tx, tz, ax, az, n)
        elif (d - e) % 3 == 0:
            d = (d - e) // 3
            (tx, tz) = point_add(ax, az, bx, bz, cx, cz, n)
            (cx, cz) = point_add(cx, cz, ax, az, bx, bz, n)
            (bx, bz, tx, tz) = (tx, tz, bx, bz)
            (tx, tz) = point_double(ax, az, n, a24)
            (ax, az) = point_add(ax, az, tx, tz, ax, az, n)
        else:
            e //= 2
            (cx, cz) = point_add(cx, cz, bx, bz, ax, az, n)
            (bx, bz) = point_double(bx, bz, n, a24)
    assert d == e == 1
    (x, z) = point_add(ax, az, bx, bz, cx, cz, n)
    return (x, z)
