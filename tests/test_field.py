from arithm.field import Field
from arithm.binary_field import BinaryField


def test_field():
    F = Field(10007)
    a = F.rand()
    b = F.rand()

    assert (a * b).val == (a.val * b.val) % F.mod
    assert a / a == F(1)


def test_binary_field():
    F = BinaryField(8, 0x11B)
    a = F.rand()

    assert a / a == F(1)
