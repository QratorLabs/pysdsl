import random
import pytest

from pysdsl import SortedIntStack


@pytest.mark.parametrize("max_value", [10, 1000, 100000])
def test_push_pop_top_empty(max_value):
    s = SortedIntStack(max_value=max_value)
    a = []

    assert s.empty()
    assert s.size == 0
    assert len(s) == 0

    # generate a sorted list of different integers less than max_value 
    a = list(range(max_value))
    random.shuffle(a)
    a = sorted(a[:max_value//2])

    for i, elem in enumerate(a):
        s.push(elem)
        assert s.top() == elem
        assert not s.empty()
        assert len(s) == i + 1
        assert s.size == i + 1

    for _ in range(max_value // 2):
        assert len(a) == len(s)
        assert len(a) == s.size
        assert not s.empty()
        assert a[-1] == s.top()
        assert a[-1] == s.pop()
        a.pop()

    assert s.empty()
    assert s.size == 0
    assert len(s) == 0


def test_copy_assign():
    s1 = SortedIntStack(max_value=1000)

    for val in (10, 100, 1000):
        s1.push(val)

    s2 = SortedIntStack(s1)

    assert len(s1) == len(s2)

    while len(s1) > 0:
        assert s1.pop() == s2.pop()


def test_errors():
    s = SortedIntStack(max_value=1000)
    
    with pytest.raises(IndexError, match="top from empty stack"):
        s.top()
    
    with pytest.raises(IndexError, match="pop from empty stack"):
        s.pop()

    with pytest.raises(ValueError, match="elements have to be pushed in strictly increasing order"):
        s.push(1)
        s.push(0)
