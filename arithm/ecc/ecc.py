from __future__ import annotations
from ..field import FieldElement
from .mults import ml
from typing import List, Union

class Point:
    """Point on a weierstrass-form elliptic curve"""

    def __init__(self, curve, x: FieldElement, y: FieldElement, z=None):
        self.x = x
        self.y = y
        if z is None:
            self.z = x.field(1)
        else:
            self.z = z
        self.curve = curve
        if curve is None:
            raise ValueError("Curve undefined")
        else:
            self.a = curve.a

    def __repr__(self):
        return f"({self.x} : {self.y} : {self.z})"

    def __eq__(self, Q: Point):
        A = self.to_affine()
        B = Q.to_affine()
        return A.x == B.x and A.y == B.y

    def __neg__(self):
        return Point(self.curve, self.x, -self.y, self.z)

    def __add__(self, Q: Point) -> Point:
        return self.j_add(Q)

    def __sub__(self, Q: Point) -> Point:
        return self.j_add(-Q)

    def __mul__(self, s: Union[FieldElement, int]) -> Point:
        if isinstance(s, FieldElement):
            x = s.val
        else:
            x = s
        if x == 2:
            return self.j_dbl()
        return ml(x, self).to_affine()

    def __rmul__(self, s: Union[FieldElement, int]):
        return self.__mul__(s)

    def complete_add_unsafe(self, other: Point):
        """Complete addition (not constant time)"""
        if self.is_at_infinity():
            return other
        if other.is_at_infinity():
            return self

        zz_s = self.z * self.z
        zz_o = other.z * other.z
        WX_s = self.x * zz_o
        WX_o = other.x * zz_s

        if WX_s == WX_o:
            ## Compare scaled y-coordinates
            WY_s = self.y * zz_o * other.z
            WY_o = other.y * zz_s * self.z
            if WY_s == WY_o:
                return 2 * self
            F = self.x.field
            return Point(self.curve, F(0), F(1), F(0))
        return self + other

    def j_dbl(self) -> Point:
        """Jacobian point doubling"""
        x, y, z = self.x, self.y, self.z

        yy = y * y
        a = x * yy
        a = a + a
        a = a + a
        b = x * x
        zz = z * z
        b = b + b + b + self.a * zz * zz

        xx = b * b - a - a
        yy = yy + yy
        yy = yy * yy
        yy = b * (a - xx) - yy - yy
        zz = y * z
        zz = zz + zz

        return Point(self.curve, xx, yy, zz)

    def j_add(self, Q: Point) -> Point:
        """Jacobian point addition"""
        x1, y1, z1 = self.x, self.y, self.z
        x2, y2, z2 = Q.x, Q.y, Q.z

        z1sq = z1 * z1
        z2sq = z2 * z2
        a = x1 * z2sq
        b = x2 * z1sq
        c = y1 * z2sq * z2
        d = y2 * z1sq * z1
        e = b - a
        f = d - c

        ee = e * e
        t = a * ee
        eee = ee * e
        xx = f * f - t - t - eee
        yy = f * (t - xx) - c * eee
        zz = z1 * z2 * e
        return Point(self.curve, xx, yy, zz)

    def dblu(self) -> Point:
        t0 = self.a
        t1 = self.x
        t2 = self.y
        t3 = t2 + t2
        t2 = t2 * t2
        t4 = t1 + t2
        t4 = t4 * t4
        t5 = t1 * t1
        t4 = t4 - t5
        t2 = t2 * t2
        t4 = t4 - t2
        t1 = t4 + t4
        t0 = t0 + t5
        t5 = t5 + t5
        t0 = t0 + t5
        t4 = t0 * t0
        t5 = t1 + t1
        t4 = t4 - t5
        t2 = t2 + t2
        t2 = t2 + t2
        t2 = t2 + t2
        t5 = t1 - t4
        t5 = t5 * t0
        t5 = t5 - t2
        return [
            Point(self.curve, t1, t2, t3),
            Point(self.curve, t4, t5, t3),
        ]

    def dblu_r(self) -> Point:
        x, y, z = self.x, self.y, self.z
        N = z * z
        E = y * y
        B = x * x
        L = E * E
        S = x.field(2) * ((x + E) ** 2 - B - L)
        M = x.field(3) * B + self.a * N**2
        xx = M * M - S - S
        zz = (y + z) ** 2 - E - N
        yy = M * (S - xx) - x.field(8) * L
        return Point(self.curve, xx, yy, zz)

    def dblu_z(self) -> List[Point]:
        A = self.a
        px, py, z = self.x, self.y, self.z
        t0 = py + z
        t0 = t0 * t0
        t1 = py * py
        t2 = z * z
        z = t0 - t1
        z = z - t2
        t0 = px + t1
        t3 = t0 * t0
        t1 = t1 * t1
        t0 = px * px
        t3 = t3 - t0
        t3 = t3 - t1
        t3 = t3 + t3
        t1 = t1 + t1
        t1 = t1 + t1
        t1 = t1 + t1
        qy = t2 * t2
        t2 = A * qy
        t2 = t2 + t0
        t2 = t2 + t0
        t2 = t2 + t0
        t0 = z * z
        qx = px
        qx = t0 * px
        qy = t0 * z
        qy = qy * py
        t0 = t2 * t2
        t0 = t0 - t3
        px = t0 - t3
        t0 = t3 - px
        py = t2 * t0
        py = py - t1

        return [Point(self.curve, qx, qy, z), Point(self.curve, px, py, z)]

    def zaddc(self, Q: Point) -> List[Point]:
        x1, y1, z = self.x, self.y, self.z
        x2, y2 = Q.x, Q.y
        c = (x1 - x2) ** 2
        w1 = x1 * c
        w2 = x2 * c
        d = (y1 - y2) ** 2
        a1 = y1 * (w1 - w2)
        x3 = d - w1 - w2
        y3 = (y1 - y2) * (w1 - x3) - a1
        z3 = z * (x1 - x2)
        d_ = (y1 + y2) ** 2
        x3_ = d_ - w1 - w2
        y3_ = (y1 + y2) * (w1 - x3_) - a1
        return [Point(self.curve, x3, y3, z3), Point(self.curve, x3_, y3_, z3)]

    def zaddu(self, Q: Point) -> List[Point]:
        x1, y1, z = self.x, self.y, self.z
        x2, y2 = Q.x, Q.y
        c = (x1 - x2) ** 2
        w1 = x1 * c
        w2 = x2 * c
        d = (y1 - y2) ** 2
        a1 = y1 * (w1 - w2)
        x3 = d - w1 - w2
        y3 = (y1 - y2) * (w1 - x3) - a1
        z3 = z * (x1 - x2)
        return [Point(self.curve, x3, y3, z3), Point(self.curve, w1, a1, z3)]

    def to_affine(self) -> Point:
        """Convert this point to affine representation (x,y,1)"""
        iz = ~self.z
        return Point(self.curve, self.x * iz**2, self.y * iz**3)

    def is_at_infinity(self) -> bool:
        """whether this point is 'zero'"""
        return self.z.val == 0
