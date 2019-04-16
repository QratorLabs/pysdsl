import array
import pysdsl


def test_intvector():
    v = pysdsl.IntVector(1024 * 1024 * 256)
    assert v.size_in_mega_bytes < 2049
    v.set_to_id()
    assert v.width == 64
    v.bit_compress()
    assert v.width == 28
    assert v.size_in_mega_bytes < 900

    v = pysdsl.IntVector([3, 2, 1, 0, 2, 1, 3, 4, 1, 1, 1, 3, 2, 3])
    v.bit_compress()
    assert sorted(v)[2] == 1
