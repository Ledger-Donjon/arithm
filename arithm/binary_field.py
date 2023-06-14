""" Field arithmetic modulo a power of two """
from secrets import randbits


class BinaryField:
    """Field modulo 2^n"""

    def __init__(self, n, mod):
        self.n = n
        self.mod = mod
        self.msb = 1 << (n - 1)

    def __repr__(self):
        return f"Binary Field of degree {self.n} and mod {self.mod:x}"

    def __call__(self, val):
        return BinaryFieldElement(val, self)

    def rand(self):
        return BinaryFieldElement(randbits(self.n), self)


class BinaryFieldElement:
    """An element belonging to a `BinaryField`"""

    def __init__(self, val, binfield):
        self.field = binfield
        if not isinstance(val, int):
            raise ValueError(f"{type(val)} is not an int")
        self.val = val

    def __repr__(self):
        return f"0x{self.val:x}"

    def __add__(self, other):
        return BinaryFieldElement(self.val ^ other.val, self.field)

    def __sub__(self, other):
        return BinaryFieldElement(self.val ^ other.val, self.field)

    def __neg__(self):
        return self

    def __eq__(self, other):
        return self.val == other.val

    def __neq__(self, other):
        return self.val != other.val

    def xtimes(self, a):
        """Helper function for multiplication by X"""
        m = self.field.msb
        return (a << 1) ^ (self.field.mod * (a & m == m))

    def __mul__(self, other):
        r = 0
        s = self.val
        if isinstance(other, BinaryFieldElement):
            t = other.val
        else:
            t = other
        for _ in range(self.field.n):
            if t & 1:
                r ^= s
            s = self.xtimes(s)
            t >>= 1
        return BinaryFieldElement(r, self.field)

    def __pow__(self, exp):
        ex = exp
        r, b = self.field(1), self
        while ex:
            if ex & 1:
                r = r * b
            b = b * b
            ex >>= 1
        return r

    def __invert__(self):
        if self.val == 0:
            raise ValueError("Trying to invert 0")
        return self ** ((1 << self.field.n) - 2)

    def __truediv__(self, other):
        return self * ~other
