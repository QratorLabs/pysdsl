import pysdsl
import pytest


@pytest.mark.parametrize("Type", pysdsl.all_wavelet_trees)
def test_wavelet(Type):
    a = Type([3, 2, 1, 0, 2, 1, 3, 4, 1, 1, 1, 3, 2, 3])
    assert a.select(2, 3) == 6 or a.select(2, 3) == 49  # 2nd is for Huffman
