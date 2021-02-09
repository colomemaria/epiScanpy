# this is to test whether episcanpy can be imported

import pytest
import episcanpy

testdata = [
    (1,2),
    (2,3),
]

def inc(x):
    return x + 1

@pytest.mark.parametrize("a,expected", testdata)
def test_simpletest(a, expected):
    assert inc(a) == expected
