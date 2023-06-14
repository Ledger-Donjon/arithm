import hashlib
from ..field import Field
from .ecc import Point
from .edwards import EdwardsPoint


class WeierstrassCurve:
    def __init__(self, p, a, b, order, gx, gy):
        F = Field(p)
        self.a = F(a)
        self.b = F(b)
        self.G = Point(self, F(gx), F(gy))
        self.order = order
        self.zero = Point(self, F(0), F(1), F(0))

    def is_on_curve(self, p):
        return p.y * p.y == p.x**3 + p.x * self.a + self.b


secp256k1 = WeierstrassCurve(
    0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F,
    0,
    7,
    0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141,
    0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798,
    0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8,
)

secp521r1 = WeierstrassCurve(
    (1 << 521) - 1,
    -3,
    0x0051953EB9618E1C9A1F929A21A0B68540EEA2DA725B99B315F3B8B489918EF109E156193951EC7E937B1652C0BD3BB1BF073573DF883D2C34F1EF451FD46B503F00,
    0x01FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFA51868783BF2F966B7FCC0148F709A5D03BB5C9B8899C47AEBB6FB71E91386409,
    0xC6858E06B70404E9CD9E3ECB662395B4429C648139053FB521F828AF606B4D3DBAA14B5E77EFE75928FE1DC127A2FFA8DE3348B3C1856A429BF97E7E31C2E5BD66,
    0x11839296A789A3BC0045C8A5FB42C7D1BD998F54449579B446817AFBD17273E662C97EE72995EF42640C550B9013FAD0761353C7086A272C24088BE94769FD16650,
)

# Edwards : x^2 + y^2 = c^2.(1 + x^2.y^2)
# Twisted : a.x^2 + y^2 = 1 + d.x^2.y^2
class TwistedEdwardsCurve:
    def __init__(self, q, order, d, a):
        self.F = Field(q)
        self.q = q
        self.r = order
        self.d = d
        self.a = self.F(a)
        self.d2 = d + d
        self.zero = EdwardsPoint(self, self.F(0), self.F(1), self.F(1))

    # https://tools.ietf.org/html/rfc8032  p.21
    # Compute corresponding x-coordinate, with low bit corresponding to
    # sign, or return None on failure
    def recover_x(self, y, sign):
        """Recover x coordinate from y coordinate and sign bit"""
        F = self.F
        modp_sqrt_m1 = F(2) ** ((self.q - 1) // 4)
        x2 = (y * y - F(1)) / (F(self.d) * y * y + F(1))
        if x2 == F(0):
            if sign:
                return None
            else:
                return 0

        # Compute square root of x2
        x = x2 ** ((self.q + 3) // 8)
        if x * x - x2 != F(0):
            x = x * modp_sqrt_m1
        if x * x - x2 != F(0):
            return None

        if (x.val & 1) != sign:
            x = -x
        return x


# Ed25519 : -x^2 + y^2 = 1 - (121665/121666)x^2.y^2
ed25519 = TwistedEdwardsCurve(
    (1 << 255) - 19,
    (1 << 252) + 27742317777372353535851937790883648493,
    0x52036CEE2B6FFE738CC740797779E89800700A4D4141D8AB75EB4DCA135978A3, # d = -F(121665) / F(121666)
    (1 << 255) - 19 - 1,
)


def secret_expand(secret):
    """Secret expansion according to rfc8032"""
    if len(secret) != 32:
        raise ValueError("Bad size of private key")
    h = hashlib.sha512(secret).digest()
    a = int.from_bytes(h[:32], "little")
    a &= (1 << 254) - 8
    a |= 1 << 254
    return (a, h[32:])


def point_decompress(s):
    """Point decompression according to rfc8032"""
    if len(s) != 32:
        raise ValueError("Invalid input length for decompression")
    y = int.from_bytes(s, "little")
    sign = y >> 255
    y &= (1 << 255) - 1
    y = ed25519.F(y)
    x = ed25519.recover_x(y, sign)
    if x is None:
        return None
    else:
        return EdwardsPoint(ed25519, x, y)

