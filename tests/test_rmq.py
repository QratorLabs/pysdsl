import random
import pytest

from pysdsl import Int8Vector, Int64Vector
from pysdsl import rmq_sada, rmq_sct




def _test_rmq(rmq_class, container_class, target):
    cont = container_class(50)
    a = []

    for i in range(50):
        a.append(random.randint(0, 50))
        cont[i] = a[i]

    rmq = rmq_class(cont)

    for _ in range(50):
        l = random.randint(0, 49)
        r = random.randint(l, 49)

        i = rmq(l, r)

        assert cont[i] == target(a[l:r+1])



@pytest.mark.parametrize('container_class', (Int8Vector, Int64Vector))
def test_rmq_sada(container_class):
    _test_rmq(rmq_sada['Min'], container_class, min)
    _test_rmq(rmq_sada['Max'], container_class, max)


@pytest.mark.parametrize('container_class', (Int8Vector, Int64Vector))
def test_rmq_sct(container_class):
    _test_rmq(rmq_sct['Min'], container_class, min)
    _test_rmq(rmq_sct['Max'], container_class, max)



@pytest.mark.parametrize("rmq_class", list(rmq_sada.values()) + list(rmq_sct.values()))
def test_rmq_on_container_remove(rmq_class):
    cont = Int64Vector(3)
    
    cont[0] = 1
    cont[1] = 2
    cont[2] = 3
    
    rmq = rmq_class(cont)

    ans1 = rmq(0, 2)

    cont[0] = 3
    cont[1] = 2
    cont[2] = 1

    ans2 = rmq(0, 2)
    assert ans1 == ans2

    cont = None

    ans2 = rmq(0, 2)
    assert ans1 == ans2
