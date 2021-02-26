# this is to test whether episcanpy can be imported

import pytest
import episcanpy
import episcanpy.api
import episcanpy.count_matrix
import episcanpy.functions
import episcanpy.plotting
import episcanpy.preprocessing
import episcanpy.tools

testdata = [
    (1,2),
    (2,3),
]

def inc(x):
    return x + 1

@pytest.mark.parametrize("a,expected", testdata)
def test_simpletest(a, expected):
    assert inc(a) == expected
