# Arithm

A small arithmetic toolbox.
Built mostly for testing out various scalar multiplication algorithms (see `./ecc/mults.py`) without having to use Sage.

Contains abstractions for Field (and binary field) arithmetic, with support for Weierstrass curves and ed25519.

## Installation

`python setup.py install` or `python setup.py develop` if you need to do modifications

## Usage

### Binary Field usage

Creating a binary field requires providing an irreducible polynomial of degree matching the size of the given binary field.

```python
from arithm.binary_field import BinaryField

BinaryField(8, 0x11b)
```

defines the AES GF(2^8) field, `0x11b` being the binary representation of X^8 + X^4 + X^3 + X + 1

```python
BinaryField(6, 0x7f)
```

defines GF(2^6) defined by X^6 + X^5 + X^4 + X^3 + X^2 + X + 1

### Field usage

```python
from arithm.field import Field

Field(10007)

a = F(57)
```