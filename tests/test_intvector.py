import pysdsl
import pytest


def test_intvector():
    v = pysdsl.IntVector(1024 * 1024 * 256)
    assert int(v.size_in_mega_bytes) == 2048
    v.set_to_id()
    assert v.width == 64
    v.bit_compress()
    assert v.width == 28
    assert v.size_in_mega_bytes < 900


@pytest.mark.parametrize("Type", [pysdsl.Int4Vector, pysdsl.Int8Vector,
                                  pysdsl.Int16Vector, pysdsl.Int24Vector,
                                  pysdsl.Int32Vector, pysdsl.Int64Vector])
def test_intNvector(Type):
    v = Type([3, 2, 1, 0, 2, 1, 3, 4, 1, 1, 1, 3, 2, 3])
    assert sorted(v)[2] == 1
    v.bit_resize(5)
    assert v.size_in_bytes == 16


def test_bitvector():
    v_len = 10
    v = pysdsl.BitVector(v_len)
    for i in range(v_len):
        v[i] = i % 2
    v.flip()
    assert v.width == 1
    assert v.max() == 1
    assert v.min() == 0
    v = pysdsl.BitVector([1, 0, 1])
    assert v.bit_size == 3


@pytest.mark.parametrize("Type", pysdsl.all_immutable_bitvectors)
def test_immutable_bitvector(Type):
    v = Type(pysdsl.BitVector([0, 1, 0, 1, 0, 1, 0, 1, 0, 1]))
    assert v.size == 10
    assert v.get_int(0, v.size) == 682
    assert v.max() == 1
    assert v.min() == 0
