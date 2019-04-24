import pysdsl
import pytest


@pytest.mark.parametrize("Type", pysdsl.all_compressed_integer_vectors)
def test_encvector(Type):
    if Type == pysdsl.pysdsl.DirectAccessibleCodesVectorDP or \
       Type == pysdsl.pysdsl.DirectAccessibleCodesVectorDPRRR:
        v = Type()
        assert v.size_in_bytes != 0
        return
    v = Type([3, 2, 1, 0, 2, 1, 3, 4, 1, 1, 1, 3, 2, 3])
    assert not v.is_sorted()
    assert v.sum() == 27
    assert v.minmax() == (0, 4)
    assert v.size_in_bytes < 200
