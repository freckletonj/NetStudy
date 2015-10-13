'''

Normalize a Connection Matrix

'''

import numpy as np
from itertools import  permutations

# CM NORMALIZATION #############################################################################

def permute_cm(cm, perm):
    '''
    permute a connection matrix according to a permutation key,
    eg the identity key: `(0,1,2)`
    or some other key such as: `(2,0,1)`
    '''
    
    num_nodes = len(perm)

    np_cm = np.array(cm)
    np_perm = np.array(perm)
    
    perm_cm = np_cm[np_perm]
    perm_cm = perm_cm.T[np_perm].T
    return perm_cm

def normalize_cm(cm):
    '''
    return the normalized connectivity matrix
    permute to all possible states, and return the lowest ordered cm
    '''
    
    num_nodes = len(cm[0])
    np_cm = np.array(cm)

    cm_perms = []
    perms = tuple(permutations(range(num_nodes)))
    for perm in perms:
        cm_perms.append(permute_cm(np_cm, perm).tolist())

    cm_perms.sort()

    return cm_perms[0]
