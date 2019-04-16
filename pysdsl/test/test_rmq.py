from pysdsl import IntVector, Int16Vector

from pysdsl import rmq_sada, rmq_sct, rmq_sparse_tables

import random


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

        assert(cont[i] == target(a[l:r+1]))


def test_rmq_sada():
    _test_rmq(rmq_sada['Min'], IntVector, min)
    _test_rmq(rmq_sada['Max'], IntVector, max)
    _test_rmq(rmq_sada['Min'], Int16Vector, min)
    _test_rmq(rmq_sada['Max'], Int16Vector, max)

def test_rmq_sct():
    _test_rmq(rmq_sct['Min'], IntVector, min)
    _test_rmq(rmq_sct['Max'], IntVector, max)
    _test_rmq(rmq_sct['Min'], Int16Vector, min)
    _test_rmq(rmq_sct['Max'], Int16Vector, max)

def test_rmq_sparse_table():
    _test_rmq(rmq_sparse_tables['Min_in_IntVector'], IntVector, min)
    _test_rmq(rmq_sparse_tables['Max_in_IntVector'], IntVector, max)
    _test_rmq(rmq_sparse_tables['Min_in_Int16Vector'], Int16Vector, min)
    _test_rmq(rmq_sparse_tables['Max_in_Int16Vector'], Int16Vector, max)
