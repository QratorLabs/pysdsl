import pysdsl
import pytest


@pytest.mark.parametrize("Type", list(pysdsl.wavelet_tree_int.values())
                         + list(pysdsl.wavelet_matrix_int.values())
                         + list(pysdsl.wavelet_tree_huffman_int.values())
                         + list(pysdsl.wavelet_tree_balanced_int.values())
                         + list(pysdsl.wavelet_tree_hu_tucker_int.values())
                         + [pysdsl.WaveletTreeGMRrankselect,
                            pysdsl.WaveletTreeGMRrankselectEnc,
                            pysdsl.WaveletTreeGolynskiMunroRao,
                            pysdsl.WaveletTreeGolynskiMunroRaoEnc,
                            pysdsl.WaveletTreeAP])
def test_wavelet(Type):
    a = Type([3, 2, 1, 0, 2, 1, 3, 4, 1, 1, 1, 3, 2, 3])
    assert a.select(2, 3) == 6


@pytest.mark.parametrize("Type", list(pysdsl.wavelet_tree_huffman.values())
                         + list(pysdsl.wavelet_tree_hu_tucker.values())
                         + list(pysdsl.wavelet_tree_balanced.values()))
def test_huffman_wavelet(Type):
    a = Type(pysdsl.BitVector([1, 0, 1, 0, 1, 0, 1, 0, 1, 0]))
    assert a.select(1, 0) == 2
