import pysdsl
import pytest


@pytest.mark.parametrize("Type", [pysdsl.SuffixArraySadakaneInt,
                                  pysdsl.SuffixArrayWaveletTreeInt])
def test_int_suffixarray(Type):
    a = Type([3, 2, 1, 5, 2, 1, 3, 4, 1, 1, 1, 3, 2, 1])
    assert a.count([3, 2, 1]) == 2
    assert a.count([2, 1]) == 3
    assert a.count([1, 2]) == 0
    assert a.count([1, 1]) == 2
    assert a.count([1]) == 6
    assert a.count([]) == 16
    assert a.text[5] == 2
    assert a.sigma == 7


@pytest.mark.parametrize("Type", [pysdsl.SuffixArrayWaveletTree,
                                  pysdsl.SuffixArraySadakane,
                                  pysdsl.SuffixArrayBitcompressed])
def test_char_suffixarray(Type):
    a = Type("abracadabra")
    assert a.count("abr") == 2
    assert a.count("a") == 5
    assert a.count("dab") == 1
    assert a.count("brac") == 1
    assert a.count("bra") == 2
    assert a.count("") == 12
    assert a.count("aba") == 0
    assert chr(a.text[5]) == "a"
    assert a.sigma == 6
