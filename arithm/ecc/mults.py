from itertools import zip_longest
from random import getrandbits


def bits(k):
    """returns a list of the bits of 'k' (LSB first)"""
    return [(k >> i) & 1 for i in range(k.bit_length())]


def bits_const(k, n):
    """returns a list of the bits of 'k' zero-padded up to 'n' bits (LSB first)"""
    return [(k >> i) & 1 for i in range(n)]


def bits_double(k, r):
    """returns a list of the bits of 'k' and 'r' (LSB first)"""

    def list2int_rev(x):
        return (x[1] << 1) | x[0]

    return list(map(list2int_rev, zip_longest(bits(k), bits(r), fillvalue=0)))


def coz_ml(k, P):
    """CoZ Montgomery Ladder"""
    R = P.dblu()
    for b in bits(k)[::-1][1:]:
        R[1 - b], R[b] = R[b].zaddc(R[1 - b])
        R[b], R[1 - b] = R[1 - b].zaddu(R[b])
    return R[0]


def ml(k, P):
    """Montgomery Ladder"""
    R = [P, 2 * P]
    for b in bits(k)[::-1][1:]:
        R[1 - b] += R[b]
        R[b] = 2 * R[b]
    return R[0]


def straus(k, P, r, Q):
    """Straus-Shamir trick for double-base multiplication"""
    s = bits_double(k, r)[::-1]
    R = [None, P, Q, P + Q]
    B = R[s[0]]
    for x in s[1:]:
        B = 2 * B
        if x != 0:
            B = B.complete_add_unsafe(R[x])
    return B


def ml_const(k, P, n):
    """Fixed-length Montgomery Ladder"""
    # always compute n-k
    e = n - k

    # do the following in constant time
    if k.bit_length() != n.bit_length():
        k_ = e
        # computing -P is 'just' computing p-y
        Q = -P
    else:
        k_ = k
        # perform a dummy subtraction somewhere so that
        # this branch matches the 'if' branch in behaviour
        Q = P

    R = [Q, 2 * Q]
    i = 0
    for b in bits_const(k_, n.bit_length())[::-1][1:]:
        R[1 - b] += R[b]
        R[b] = 2 * R[b]
        i += 1
    # checking that we always perform
    # the same number of iterations
    assert i == n.bit_length() - 1
    return R[0]


def r2l_daa(k, n, P):
    """Right-to-left double-and-add"""
    # Avoiding adding Q to itself in the loop :
    # - Initialize R0, R1, with P
    # - skip first bit so that the accumulator B is necessarily an even multiple of P
    # - Rx is always an odd multiple, therefore always different from B.
    R = [P, P]  # [ Fake, Result ]
    B = 2 * P  # We skip the first bit
    kbits = bits_const(k, n)
    for b in kbits[1:]:
        R[b] += B
        B = 2 * B

    b = 1 - kbits[0]  # If first bit was 0, subtract P
    R[b] = R[b] - P  # If first bit was set, then result is correct (don't subtract P)
    return R[1]


def r2l_daa_w(k, n, P, w):
    """Windowed Right-to-left double-and-add"""
    R = [P for i in range(1 << w)]
    B = 2 * P
    kbits = bits_const(k, n + 1 + (w - n % w))
    for i in range(1, len(kbits), w):
        b = 0
        for j in range(w - 1, -1, -1):
            b <<= 1
            b |= kbits[i + j]
        R[b] += B
        for j in range(w):
            B = 2 * B

    for i in range(2, 1 << w):
        R[1] += i * R[i]

    R[1] -= 27 * P
    b = 1 - kbits[0]
    R[b] = R[b] - P
    return R[1]


def r2l_daa_point_blinding(k, n, P):
    """Right-to-left double-and-add with point blinding"""
    # Assuming we have precomputed some random point that
    # is different from any intermediate value and avoids adding
    # two identical points during the algorithm.
    # Picking a point at random has no guarantee beyond statistical.

    # Below can be replaced with a much cheaper point mapping algorithm
    # like SWU or Icart
    r = getrandbits(32)
    Q = ml(r, P)

    R = [Q, Q]  # [ Fake, Result ]
    B = P
    for b in bits_const(k, n):
        R[b] += B
        B = 2 * B
    return R[1] - Q
