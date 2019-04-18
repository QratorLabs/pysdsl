from pysdsl import *


def test_encvector():
    for T in [EncVectorEliasDelta, EncVectorEliasGamma, EncVectorFibonacci,
              EncVectorComma4, EncVectorComma2]:
        v = T([3, 2, 1, 0, 2, 1, 3, 4, 1, 1, 1, 3, 2, 3])
        assert not v.is_sorted()
        assert v.sum() == 27
        assert v.minmax() == (0, 4)
        assert v.size_in_bytes < 200
