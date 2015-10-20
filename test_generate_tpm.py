import numpy as np
from generate_tpm import generate_tpm, OR, AND, XOR

np.testing.assert_almost_equal( generate_tpm([OR, AND, XOR],  [[0,1,1], [1,0,1],[1,1,0]]),
                                    [[0,0,0],[0,0,1],[1,0,1],[1,0,0],[1,0,0],[1,1,1],[1,0,1],[1,1,0]])
