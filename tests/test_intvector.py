from pysdsl import *


def test_intvector():
    v = IntVector(1024 * 1024 * 256)
    assert v.size_in_mega_bytes < 2049
    v.set_to_id()
    assert v.width == 64
    v.bit_compress()
    assert v.width == 28
    assert v.size_in_mega_bytes < 900

    for T in [IntVector, Int4Vector, Int8Vector, Int16Vector, Int24Vector,
              Int32Vector, Int64Vector]:
        v = IntVector([3, 2, 1, 0, 2, 1, 3, 4, 1, 1, 1, 3, 2, 3])
        v.bit_compress()
        assert sorted(v)[2] == 1

    v_len = 10
    bit_v = BitVector(v_len)
    for i in range(v_len):
        bit_v[i] = i % 2
    bit_v.flip()
    assert bit_v.width == 1
    assert bit_v.max() == True
    assert bit_v.min() == False
    bit_v = BitVector([1, 0, 1])
    assert bit_v.bit_size == 3
