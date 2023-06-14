""" Field arithmetic modulo a prime number"""
from secrets import randbits
from sympy.ntheory.primetest import isprime


def exgcd(a, b):
    ## Handbook of Elliptic and Hyperelliptic Curve Cryptography 10.6.1
    """Return (u,v,d) such that a.u + b.v = d = gcd(a,b)"""
    aa, bb = b, a
    ua, ub = 0, 1
    while bb != 0:
        q = aa // bb
        aa, bb = bb, aa - q * bb
        ua, ub = ub, ua - q * ub
    v = (aa - a * ua) // b
    return ua, v, aa


def invmod(x, m):
    """Helper function for inversion of `x` modulo `m`"""
    u, _, d = exgcd(x, m)
    if d != 1:
        raise ValueError(f"{x} is not invertible modulo {m}")
    return u % m


class Field:
    """Field modulo a prime number"""

    def __init__(self, mod):
        if not isprime(mod):
            print(f"Warning: {mod} does not appear to be prime")
        self.mod = mod

    def __repr__(self):
        return f"Field modulo {self.mod}"

    def __call__(self, val):
        return FieldElement(val, self)

    def rand(self):
        """Get a random element"""
        return FieldElement(randbits(self.mod.bit_length() + 64), self)


class FieldElement:
    """An element belonging to a `Field`"""

    def __init__(self, val, field):
        self.field = field
        if not isinstance(val, int):
            raise ValueError(f"{type(val)} is not supported")
        self.val = val % field.mod

    def __repr__(self):
        return hex(self.val)

    def __add__(self, other):
        assert self.field.mod == other.field.mod
        return FieldElement((self.val + other.val) % self.field.mod, self.field)

    def __sub__(self, other):
        assert self.field.mod == other.field.mod
        return FieldElement((self.val - other.val) % self.field.mod, self.field)

    def __neg__(self):
        return FieldElement(self.field.mod - self.val, self.field)

    def __eq__(self, other):
        assert self.field.mod == other.field.mod
        return self.val == other.val

    def __neq__(self, other):
        assert self.field.mod == other.field.mod
        return self.val != other.val

    def __mul__(self, other):
        if isinstance(other, FieldElement):
            assert self.field.mod == other.field.mod
            t = other.val
        elif isinstance(other, int):
            t = other
        else:
            return NotImplemented
        return FieldElement((self.val * t) % self.field.mod, self.field)

    def __invert__(self):
        return FieldElement(invmod(self.val, self.field.mod), self.field)

    def __truediv__(self, other):
        return self * invmod(other.val, self.field.mod)

    def __pow__(self, exp):
        return FieldElement(pow(self.val, exp, self.field.mod), self.field)

    def legendre(self):
        """Compute the legendre symbol"""
        k = 1
        a = self.val
        p = self.field.mod
        while p != 1:
            if a == 0:
                return 0
            v = 0
            while a % 2 == 0:
                v += 1
                a //= 2
            pmod4 = p % 4
            pmod8 = p % 8
            if v % 2 == 1 and (pmod8 == 3 or pmod8 == 5):
                k = -k
            if a % 4 == 3 and pmod4 == 3:
                k = -k
            r = a
            a = p % r
            p = r
        return k

    def sqrt(self):
        """Compute the square root"""
        ## TODO: does not do Tonnelli-Shanks yet
        if self.legendre() == 1 and self.field.mod % 4 == 3:
            return self ** ((self.field.mod + 1) // 4)
        else:
            raise ValueError(f"{self} is not a square")
