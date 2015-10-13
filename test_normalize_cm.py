import unittest
import normalize_cm
from numpy.testing import assert_array_equal

def test_permute_cm():
    cm = [[1,2,3],
          [4,5,6],
          [7,8,9]]

    cm201 = [[9,7,8],
             [3,1,2],
             [6,4,5],]

    assert_array_equal(normalize_cm.permute_cm(cm, (2,0,1)), cm201)
