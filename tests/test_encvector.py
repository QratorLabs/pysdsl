import pysdsl
import pytest


@pytest.mark.parametrize("Type", pysdsl.enc_vector.values()
                         + pysdsl.variable_length_codes_vector.values()
                         + [pysdsl.DirectAccessibleCodesVector4,
                            pysdsl.DirectAccessibleCodesVector8,
                            pysdsl.DirectAccessibleCodesVector16,
                            pysdsl.DirectAccessibleCodesVector63])
def test_encvector(Type):
    v = Type([3, 2, 1, 0, 2, 1, 3, 4, 1, 1, 1, 3, 2, 3])
    assert not v.is_sorted()
    assert v.sum() == 27
    assert v.minmax() == (0, 4)
    assert v.size_in_bytes < 200


@pytest.mark.skip(reason="Issue #16")
@pytest.mark.parametrize("Type", [pysdsl.DirectAccessibleCodesVectorDP,
                                  pysdsl.DirectAccessibleCodesVectorDPRRR])
def test_encvectordp(Type):
    v = Type([3, 2, 1, 0, 2, 1, 3, 4, 1, 1, 1, 3, 2, 3])
    assert not v.is_sorted()
    assert v.sum() == 27
    assert v.minmax() == (0, 4)
    assert v.size_in_bytes < 200
