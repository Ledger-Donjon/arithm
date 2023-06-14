from binascii import unhexlify
from random import getrandbits
from arithm.field import Field
from arithm.ecc.edwards import EdwardsPoint
from arithm.ecc.mults import *
from arithm.ecc.curves import secp256k1, ed25519, secret_expand, point_decompress


def test_secp256k1():
    """random testing of Montgomery Ladder, CoZ ML, Right-to-left double add-always, it's variant with point blinding, and k/n-k ML"""
    P = secp256k1.G
    n = secp256k1.order

    # test r2l, const_ml's correctness vs ml
    for _ in range(10):
        k = getrandbits(256)

        ref = ml(k, P).to_affine()
        coz = coz_ml(k, P).to_affine()
        r2l = r2l_daa(k, 256, P).to_affine()
        r2lw = r2l_daa_w(k, 256, P, 3).to_affine()
        r2lb = r2l_daa_point_blinding(k, 256, P).to_affine()
        mlc = ml_const(k, P, n).to_affine()

        assert ref == coz
        assert ref == r2l
        assert ref == r2lw
        assert ref == r2lb
        assert ref == mlc


def test_straus_secp256k1():
    """Straus' trick for secp256k1"""
    P = secp256k1.G

    # test straus's correctness
    for _ in range(10):
        k = getrandbits(256)
        r = getrandbits(256)
        rq = getrandbits(32)
        Q = rq * P

        ref = k * P + r * Q
        ref = ref.to_affine()
        sts = straus(k, P, r, Q).to_affine()

        assert ref == sts


def test_straus_failures():
    """edge cases for Straus. This test asserts the point is at infinity at the end, which is actually a failure"""
    P = secp256k1.G
    n = secp256k1.order

    # 1. 2*(k[:i].P + r[:i].Q) = P, k[i] = 1, r[i] = 0
    # 2. 2*(k[:i].P + r[:i].Q) = Q, k[i] = 0, r[i] = 1
    # 3. 2*(k[:i].P + r[:i].Q) = P+Q, k[i] = 1, r[i] = 1
    # 4. all of the above but with inverse points

    # 1) k = n, r = 0 ...
    #    Q = (1/2-k[:i]).r[:i]^-1.P
    # 2) symmetrical
    # 3) Q = (1/2-k[:i]).(1/2-r[:i])^-1.P

    # 1)
    M = 31  ## arbitrary step at which the algorithm goes into infinity
    k = getrandbits(M)
    r = getrandbits(M)
    Fn = Field(n)
    e = ~Fn(r) * Fn((n + 1) // 2 - k)
    Q = e.val * P
    r <<= 1
    k <<= 1
    k += 1
    ## rest can be anything
    r <<= 256 - M - 1
    r += getrandbits(256 - M - 1)
    k <<= 256 - M - 1
    k += getrandbits(256 - M - 1)

    assert not straus(k, P, r, Q).is_at_infinity()


def test_edwards_basic():
    """ed25519 point decompression test"""
    # https://tools.ietf.org/html/rfc8032  7.1
    s = unhexlify("9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60")
    s = secret_expand(s)[0]
    sP = unhexlify("d75a980182b10ab7d54bfed3c964073a0ee172f3daa62325af021a68f707511a")

    sP = point_decompress(sP)

    # base EdwardsPoint :
    F = ed25519.F
    P = EdwardsPoint(ed25519, x=F(0), y=F(4) / F(5))
    P.x = ed25519.recover_x(P.y, 0)

    assert s * P == sP
