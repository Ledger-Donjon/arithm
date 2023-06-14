from ..field import FieldElement


class EdwardsPoint:
    def __init__(
        self,
        curve,
        x: FieldElement,
        y: FieldElement,
        z: FieldElement = None,
        t: FieldElement = None,
    ):
        self.curve = curve
        self.x = x
        self.y = y
        if z is None:
            self.z = x.field(1)
        else:
            self.z = z
        if t is None:
            self.t = self.x * self.y
        else:
            self.t = t

    def __repr__(self):
        return f"({self.x} : {self.y} : {self.z})"

    def __eq__(self, Q):
        A = self.to_affine()
        B = Q.to_affine()
        return A.x == B.x and A.y == B.y

    def __neg__(self):
        return EdwardsPoint(self.curve, -self.x, self.y, self.z, -self.t)

    def __add__(self, Q):
        return self.add(Q).to_affine()

    ## FIXME
    # def __rmul__(self, s):
    #     # r = Edwardspoint(self.curve,)
    #     r = self
    #     for i in bin(s)[2:]:
    #         # r = r.add(r)
    #         r = r.idbl()
    #         if int(i,2) == 1:
    #             r = r.add(self)
    #     return r.to_affine()

    def __rmul__(self, k):
        r = self.curve.zero
        s = self
        if isinstance(k, FieldElement):
            k_ = k.val
        else:
            k_ = k
        for i in bin(k_)[2:][::-1]:
            if int(i, 2):
                r = r.add(s)
            # t = s
            # s = t.add(s) ## FIXME
            s = s.idbl()
        return r.to_affine()

    # Using extended coordinates
    def to_affine(self):
        """Convert this point to affine representation (x,y,1)"""
        iz = ~self.z
        return EdwardsPoint(self.curve, self.x * iz, self.y * iz)

    # http://hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html#addition-add-2008-hwcd-3
    # assumes a = -1
    # unified and complete
    def add(self, Q):
        """Unified addition"""
        X1, Y1, Z1, T1 = self.x, self.y, self.z, self.t
        X2, Y2, Z2, T2 = Q.x, Q.y, Q.z, Q.t
        A = (Y1 - X1) * (Y2 - X2)
        B = (Y1 + X1) * (Y2 + X2)
        C = T1 * self.curve.d2 * T2
        D = Z1 * self.curve.F(2) * Z2
        E = B - A
        FF = D - C
        G = D + C
        H = B + A
        X3 = E * FF
        Y3 = G * H
        T3 = E * H
        Z3 = FF * G
        return EdwardsPoint(self.curve, X3, Y3, Z3, T3)

    def idbl(self):
        """Initial doubling"""
        X1, Y1, Z1 = self.x, self.y, self.z
        A = X1**2
        B = Y1**2
        C = self.curve.F(2) * Z1**2
        D = self.curve.a * A
        E = (X1 + Y1) ** 2 - A - B
        G = D + B
        FF = G - C
        H = D - B
        X3 = E * FF
        Y3 = G * H
        T3 = E * H
        Z3 = FF * G
        return EdwardsPoint(self.curve, X3, Y3, Z3, T3)
