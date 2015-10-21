import numpy as np
from generate_tpm import generate_tpm, OR, AND, XOR, GREATER_THAN


def test_generate_tpm():
    np.testing.assert_almost_equal( generate_tpm([OR, AND, XOR],  [[0,1,1], [1,0,1],[1,1,0]]),
                                    [[0,0,0],[0,0,1],[1,0,1],[1,0,0],[1,0,0],[1,1,1],[1,0,1],[1,1,0]])

def test_GREATER_THAN():
    G = GREATER_THAN(4)
    assert G([1,1,1,1,1])
    assert G([1,1,1,1])
    

    assert not G([1,1,1])
