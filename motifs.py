'''

For Studying : http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0020369
And comparing their results against the theory of Integrated Information

TODO: This code is built for deterministic networks, make it work for indeterminate ones as well

'''

import numpy as np
from itertools import combinations, permutations, product
from functools import reduce
from statistics import mean

from normalize_cm import normalize_cm
from normalize_tpm import normalize_tpm
from generate_tpm import generate_tpm, nodes_to_short, AND, OR, XOR, NULL

from IPython.core.debugger import Tracer

import pyphi

  
# GENERATE TEST CMs #############################################################################
# number of unique_cms:
# 3 nodes = 104 : calculates fast
# 4 nodes = 3044 :calculates in <1 min
# 5 nodes = scares me...


def main():
    NUM_NODES = 3

    cms = product(*[[0,1]]*(NUM_NODES**2)) # Warning: mind-numbingly inefficient (brute force), don't use a lot of nodes (like, more than 4 sounds scary but I haven't tried it)

    unique_cms = set() #there will be 104 unique cms at 3 nodes with recurrent loops

    all_activations =  [p for p in product([0,1], repeat=NUM_NODES)] # all possible activation states
    past_state = [0 for _ in range(NUM_NODES)]

    # look at all cms, add uniques to unique_cms
    for cm in cms:
        np_cm = np.array(cm)
        np_cm.shape = (NUM_NODES, NUM_NODES)
        unique_cms.add(tuple([tuple(x) for x in normalize_cm(np_cm)]))








        
    ########################################################################################
    # Condition (1)  -  No self loops (that is each node has at most 2 inputs and 2 outputs):
    # number of cms in condition 1:
    # 3 nodes = 13
    # 4 nodes = 202

    has_self_loops = lambda cm: bool(len([1 for i in range(NUM_NODES) if cm[i][i] == 1]))
    has_unconnected_node = lambda cm: bool(len(
        [1 for i in range(NUM_NODES) if
         sum(cm[i])==0 and
         sum(np.array(cm).T.tolist()[0])==0
     ]))

    cms1 = [cm for cm in unique_cms
            if not has_self_loops(cm)
            and not has_unconnected_node(cm)]

    cms1.sort()


    # 1.a) only linear threshold mechanisms: two possibilities Threshold >=1 (Copy/OR), Threshold >= 2 (AND)
    #       maybe only those nodes that actually have 2 inputs should be ANDs (not sure about this one)

    seen_tpms_by_condition = {}


    MECHANISMS = [OR, AND, NULL, ]

    # arr[x] = [mechs] where x is the sum of inputs, and [mechs] is the set of acceptable mechs
    INPUTS_SUM_TO_MECH_MAP = [[NULL],
                              [OR],
                              [OR, AND, ],]

    duplicate_tpms = 0
    
    for cm in cms1:
        #print("CONNECTIVITY MATRIX: ", cm)

        seen_tpms_by_condition[cm] = {}
        seen_tpms_by_condition[cm]['seen'] = set()
        seen_tpms_by_condition[cm]['num_concepts'] = []
        seen_tpms_by_condition[cm]['phi_concepts'] = []
        seen_tpms_by_condition[cm]['phi_network'] = []

        for nodes in product(MECHANISMS, repeat=NUM_NODES):
            #print("NODES: ", nodes_to_short(nodes))

            # ignore systems that don't match the INPUTS_SUM_TO_MECH_MAP
            ignore = False
            for sm, mechs in enumerate(INPUTS_SUM_TO_MECH_MAP):
                if not all_nodes_where_incoming_sum_is_x_are_mech_y(cm, nodes, sm, mechs):
                    ignore = True
                    break

            if ignore:
                continue

            tpm = generate_tpm(nodes, cm)
            n_tpm = normalize_tpm(tpm)
            
            # continue next loop if this tpm has already been seen on this condition
            if to_2d_tuple(n_tpm) in seen_tpms_by_condition[to_2d_tuple(cm)]['seen']:                    
                duplicate_tpms += 1
                continue
                
            seen_tpms_by_condition[cm]['seen'].add(to_2d_tuple(n_tpm))
                
            for current_state in all_activations:
                #print("CURRENT STATE: ", current_state)

                network = pyphi.Network(tpm, current_state, past_state, connectivity_matrix=cm)

                subsystem = pyphi.Subsystem(range(network.size), network)

                constellations = pyphi.compute.constellation(subsystem)
                
                seen_tpms_by_condition[cm]['num_concepts'].append(len(constellations))

                if len(constellations) > 0:
                    seen_tpms_by_condition[cm]['phi_concepts'].append(mean([x.phi for x in constellations]))

                    
                seen_tpms_by_condition[cm]['phi_network'].append(pyphi.compute.big_phi(subsystem))

                
                #if cm == cms1[12]:
                #    Tracer()()
                


                
    print("\n"*3)

    #from collections import OrderedDict
    #seen = OrderedDict(sorted(seen_tpms_by_condition.items()))

    # motifs from the paper "Motifs in the Brain", same ordering
    motifs_mapping = [
        [[0,0,0], [1,0,0], [1,0,0]],
        [[0,0,0], [1,0,0], [0,1,0]],
        [[0,0,0], [1,0,1], [0,0,0]],
        [[0,1,0], [1,0,0], [1,0,0]],
        [[0,0,0], [1,0,0], [1,1,0]],
        [[0,1,1], [1,0,0], [0,0,0]],
        [[0,0,1], [1,0,0], [0,1,0]],
        [[0,1,0], [1,0,0], [1,1,0]],
        [[0,1,1], [1,0,0], [1,0,0]],
        [[0,1,0], [1,0,1], [1,0,0]],
        [[0,0,0], [1,0,1], [1,1,0]],
        [[0,1,1], [1,0,0], [1,1,0]],
        [[0,1,1], [1,0,1], [1,1,0]],
        ]

    seen = seen_tpms_by_condition
    
    for i, m in enumerate(motifs_mapping):
        row = seen[to_2d_tuple(normalize_cm(m))]
        print('%2d'%(i+1), '\t',
              m, '\t',
              len( row['seen']), '\t',
              '%.2f' % (mean(row['num_concepts'])) if len(row['num_concepts'])>0 else 'null' , '\t',
              '%.2f' % (mean(row['phi_concepts'])) if len(row['phi_concepts'])>0 else 'null' , '\t',
              '%.2f' % (mean(row['phi_network'])) if len(row['phi_network'])>0 else 'null'

              
        )
    
    print("Duplicate TPMS:", duplicate_tpms)
    





    


    # 1.b) same as (1a) with the possibility of XORs for those nodes that
    #       actually have 2 inputs












    

    ########################################################################################
    # Condition (2)  -   All the nodes have self loops

    # 2.a) only linear threshold mechanisms: >= 1, >= 2, >= 3


    # 2.b) same as (2a) but with the possibility of XORs






#########################################################################################
# Utility Functions
#########################################################################################

def to_2d_tuple(arr):
    ''' takes a (mutable) 2D array and converts it to an (immutable) tuple '''
    return tuple(tuple(row) for row in arr)

    
def count_outgoing(cm):
    '''
    TODO: correct counts for non-deterministic connections
    TODO: IE strengths aren't 0 or 1
    
    return the count of outgoing connections for each node

    args:
      cm: a connectivity matrix / 2d array
          cm[0][x] = outgoing connections
          cm[x][0] = incoming connections

    returns:
      a list of how many outgoing connections each node has, eg:
      [2, 2, 0, 1]
    
      where node0 has 2 outgoing connections, node2 has none, etc.
    '''
    counts = [sum(row) for row in cm]
    return counts

    
def count_incoming(cm):
    '''
    TODO: correct counts for non-deterministic connections
    TODO: IE strengths aren't 0 or 1
    
    return the count of incoming connections for each node

    args:
      cm: a connectivity matrix / 2d array
          cm[0][x] = outgoing connections
          cm[x][0] = incoming connections

    returns:
      a list of how many incoming connections each node has, eg:
      [2, 2, 0, 1]
    
      where node0 has 2 incoming connections, node2 has none, etc.
    '''
    counts = reduce(lambda xs, ys: [x+y for x, y in zip(xs, ys)], cm)
    return counts

def all_nodes_where_incoming_sum_is_x_are_mech_y(cm, nodes, sm, mechs):
    # if a node has one input, it must be a COPY (IE OR)
    
    in_count = count_incoming(cm)
    for count, node in zip(in_count, nodes):
        if (count == sm):
            if (node in mechs):
                continue
            return False
    return True
    




if __name__ == "__main__":
    main()
