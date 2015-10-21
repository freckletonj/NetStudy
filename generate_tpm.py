'''

Generate a tpm given node mechanisms, and a connectivity matrix

'''



import pyphi
from itertools import product
import numpy as np

# GATES #############################################################################
AND =  lambda xs: sum(xs) == len(xs)
NAND =  lambda xs: not AND(xs)
OR =  lambda xs: sum(xs) >= 1
NOR =  lambda xs: not OR(xs)
XOR =  lambda xs: sum(xs) % 2 == 1
RANDOM =  None
MAJORITY =  lambda xs: sum(xs) > len(xs)/2
MINORITY =  lambda xs: sum(xs) <= len(xs)/2
PARITY =  lambda xs: sum(xs) % 2 == 0 #even parity
GREATER_THAN = lambda threshold: lambda xs: sum(xs) >= threshold
LESS_THAN = lambda threshold: lambda xs: sum(xs) < threshold
NULL = lambda _: False


    
def generate_tpm(nodes, cm):
    '''
    State-by-node representation
    
    Structure:
    t0: 000, 100, 010, 110, 001, ...
    t1: according to mechanisms
    '''
    starting_states = [list(reversed(x)) for x in product(*tuple([[0,1]]*len(nodes)))]
    tpm = np.zeros((2**len(nodes), len(nodes)))

    for i, state in enumerate(starting_states):
        tpm[i] = [node([state[x]*cm[x][node_num] for x in range(len(state)) if cm[x][node_num]>0])
                  for node_num, node in enumerate(nodes)]
        
        
    return tpm



    

    
def nodes_to_name(nodes):
    '''
    give it a node mechanism and it returns the name the of the node, useful for printing
    '''
    mapping = {
        AND: "AND",
        NAND: "NAND",
        OR: "OR",
        NOR: "NOR",
        XOR: "XOR",
        RANDOM: "RANDOM",
        MAJORITY: "MAJORITY",
        MINORITY: "MINORITY",
        PARITY: "PARITY",
        GREATER_THAN: "GREATER_THAN",
        LESS_THAN: "LESS_THAN",
        NULL: "NULL",
    }
        
    output = []
    for node in nodes:
        output.append(mapping[node])
    return output






    
def nodes_to_short(nodes):
    '''
    give it a node mechanism and it returns the name the of the node, useful for printing
    '''
    mapping = {
        AND: "A",
        NAND: "N",
        OR: "O",
        NOR: "Q",
        XOR: "X",
        RANDOM: "R",
        MAJORITY: "J",
        MINORITY: "M",
        PARITY: "P",
        GREATER_THAN: "G",
        LESS_THAN: "L",
        NULL: "U",
    }
        
    output = []
    for node in nodes:
        output.append(mapping[node])
    return output
