from pysdsl import *


def test_wavelet():
    for T in all_wavelet_trees:
        a = T([3, 2, 1, 0, 2, 1, 3, 4, 1, 1, 1, 3, 2, 3])
        a.size_in_bytes
