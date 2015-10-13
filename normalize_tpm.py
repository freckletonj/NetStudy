'''

Normalize a Transition Probability Matrix
    equivalent tpms normalize to the same tpm

note:
the current code runs in at least factorial time

optimizations will be added down the road

'''


import numpy as np
from itertools import permutations
from pyphi.convert import state2loli_index, loli_index2state

def input_state_permutation(num_nodes, perm):
    """
    Permute input state order according to permutation of nodes given
    """
    # in the permuted tpm the inputs nodes are ordered as in the permutation i.e. 1 2 0 
    input_states = np.array([loli_index2state(s, num_nodes) for s in range(2**num_nodes)])

    # now I have to find the permutation that would convert perm back to 0 1 2
    perm_back = [i[0] for i in sorted(enumerate(perm), key=lambda x:x[1])]
    # now permute columns of input_states using perm_back
    # used transposed cause I haven't figured out how to index columns
    permuted_input_states = input_states.T[np.array(perm_back)].T

    permuted_state_indices = [state2loli_index(permuted_input_states[s]) for s in range(2**num_nodes)]

    return permuted_state_indices


def permute_tpm(tpm, perm):
    """
    Permute state by node tpm according to permutation of nodes given
    """
    num_nodes = len(perm)
    
    perm_tpm_temp = np.zeros((2**num_nodes, num_nodes)).astype(int)
    perm_tpm = np.zeros((2**num_nodes, num_nodes)).astype(int)

    # -------- permute rows of tpm ---------
    perm_indices = input_state_permutation(num_nodes, perm)
    perm_tpm_temp = tpm[perm_indices]

    # -------- permute columns of tpm ------
    perm_tpm = perm_tpm_temp.T[np.array(perm)].T

    #import pdb; pdb.set_trace()
    return perm_tpm


def normalize_tpm(tpm):
    """Render a state by node tpm into a normal form.
    Two tpms will have the same normal form if they are
    equivalent up to a permutation of the ordering of their nodes.
    """
    # Cast to a numpy array
    tpm = np.array(tpm)
    # Get the number of nodes (length of columns)
    num_nodes = tpm.shape[1]
    # Get the permutations of the node indices
    perms = tuple(permutations(range(num_nodes)))

    tpm_perms = []
    for perm in perms:
        ptpm = permute_tpm(tpm, perm)
        tpm_perms.append(ptpm.tolist())

    tpm_perms.sort()
    return tpm_perms[0]
