import pysdsl
import pytest


@pytest.mark.parametrize("Type", [# pysdsl.SuffixArraySadakaneInt,
                                  # pysdsl.SuffixArrayWaveletTreeInt,
                                  pysdsl.SuffixArrayWaveletTree,
                                  pysdsl.SuffixArraySadakane,
                                  pysdsl.SuffixArrayBitcompressed])
def test_suffixarray(Type):
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
