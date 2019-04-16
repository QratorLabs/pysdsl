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

    v = pysdsl.Int64Vector([1, 2, 3])
    u = array.array('Q', v)
    assert u[1][2] == 2
