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
from generate_tpm import GREATER_THAN as G

from IPython.core.debugger import Tracer

import pyphi

  
# GENERATE TEST CMs #############################################################################
# number of unique_cms:
# 3 nodes = 104 : calculates fast
# 4 nodes = 3044 :calculates in <1 min
# 5 nodes = scares me...


def main():
    NUM_NODES = 3

    unique_cms = gen_unique_cms(NUM_NODES)
    all_activations =  [p for p in product([0,1], repeat=NUM_NODES)] # all possible activation states
    past_state = [0 for _ in range(NUM_NODES)]

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

    # all nodes have self loop
    motifs_plus_self_mapping = [
        [[1,0,0], [1,1,0], [1,0,1]],
        [[1,0,0], [1,1,0], [0,1,1]],
        [[1,0,0], [1,1,1], [0,0,1]],
        [[1,1,0], [1,1,0], [1,0,1]],
        [[1,0,0], [1,1,0], [1,1,1]],
        [[1,1,1], [1,1,0], [0,0,1]],
        [[1,0,1], [1,1,0], [0,1,1]],
        [[1,1,0], [1,1,0], [1,1,1]],
        [[1,1,1], [1,1,0], [1,0,1]],
        [[1,1,0], [1,1,1], [1,0,1]],
        [[1,0,0], [1,1,1], [1,1,1]],
        [[1,1,1], [1,1,0], [1,1,1]],
        [[1,1,1], [1,1,1], [1,1,1]],
    ]




        
    ########################################################################################
    # Condition (1)  -  No self loops (that is each node has at most 2 inputs and 2 outputs):
    # number of cms in condition 1:
    # 3 nodes = 13
    # 4 nodes = 202

    # returns True if there is a self-loop on any nodes
    has_self_loops = lambda cm: bool(len([1 for i in range(NUM_NODES) if cm[i][i] == 1]))

    # returns True if any node doesn't have a connection with at least one other node
    # note, this doesn't guarantee a weakly connected graph (ex A->B C->D, AB isn't connected to CD)
    has_unconnected_node = lambda cm: bool(len(
        [1 for i in range(NUM_NODES) if
         sum(cm[i])==0 and
         sum(np.array(cm).T.tolist()[i])==0
     ]))

    # all cms for condition 1
    cms1 = [cm for cm in unique_cms
            if not has_self_loops(cm)
            and not has_unconnected_node(cm)]

    # same as cms1, but each node has a self loop
    cms2 = [to_2d_list(cm) for cm in cms1]
    for i, cm in enumerate(cms2):
        for j in range(NUM_NODES):
            cm[j][j] = 1
        cms2[i] = to_2d_tuple(cm)


   

    # cms1.sort()
    # for cm in cms1:
    #     print(cm)
    # print(len(cms1))
    

    # 1.a) only linear threshold mechanisms: two possibilities Threshold >=1 (Copy/OR), Threshold >= 2 (AND)
    #       maybe only those nodes that actually have 2 inputs should be ANDs (not sure about this one)


    condition_names = ['1A', '1B', '2A', '2B']

    # Greater-Than mechanisms
    G1 = G(1)
    G2 = G(2)
    G3 = G(3)

    # sum_to_mech_map = if sum==<index of this array>, mech must be in <value @ index>
    conditions = {
        '1A':{
            'cms': cms1,
            'mechanisms': [OR, AND, NULL],
            'inputs_sum_to_mech_map': [[NULL],
                                       [OR],
                                       [OR, AND]]
        },
        '1B':{
            'cms': cms1,
            'mechanisms': [OR, AND, NULL, XOR],
            'inputs_sum_to_mech_map': [[NULL],
                                       [OR],
                                       [OR, AND, XOR]]
        },
        '2A':{
            'cms': cms2,
            'mechanisms': [G1, G2, G3, NULL],
            'inputs_sum_to_mech_map': [[NULL],
                                       [G1],
                                       [G1, G2],
                                       [G1, G2, G3]]
        },
        '2B':{
            'cms': cms2,
            'mechanisms': [G1, G2, G3, NULL, XOR],
            'inputs_sum_to_mech_map':  [[NULL],
                                        [G1],
                                        [G1, G2, XOR],
                                        [G1, G2, G3, XOR]]
        },
    }

    

    MECHANISMS = [OR, AND, NULL, ]

    # arr[x] = [mechs] where x is the sum of inputs, and [mechs] is the set of acceptable mechs
    INPUTS_SUM_TO_MECH_MAP = [[NULL],
                              [OR],
                              [OR, AND,  ],]



    for condition_name in condition_names:
        condition = conditions[condition_name]
        MECHANISMS = condition['mechanisms']
        INPUTS_SUM_TO_MECH_MAP = condition['inputs_sum_to_mech_map']
        CMS = condition['cms']
        
        print('\n\nCondition: ', condition_name)
        
        seen_tpms_by_condition = {}
        for cm in CMS:
            #print("CONNECTIVITY MATRIX: ", cm)
            #continue

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
                    continue

                seen_tpms_by_condition[cm]['seen'].add(to_2d_tuple(n_tpm))

                for current_state in all_activations:
                    #print("CURRENT STATE: ", current_state)

                    network = pyphi.Network(tpm, current_state, past_state, connectivity_matrix=cm)
                    subsystem = pyphi.Subsystem(range(network.size), network)
                    constellations = pyphi.compute.constellation(subsystem)


                    # Record stats of interest
                    seen_tpms_by_condition[cm]['num_concepts'].append(len(constellations))
                    if len(constellations) > 0:
                        seen_tpms_by_condition[cm]['phi_concepts'].append(mean([x.phi for x in constellations]))
                    else:
                        seen_tpms_by_condition[cm]['phi_concepts'].append(0)
                    seen_tpms_by_condition[cm]['phi_network'].append(pyphi.compute.big_phi(subsystem))


                    #if cm == cms1[12]:
                    #    Tracer()()

    
        #from collections import OrderedDict
        #seen = OrderedDict(sorted(seen_tpms_by_condition.items()))


        seen = seen_tpms_by_condition

        printer = lambda i, m, row: ('\'%2d'%(i+1), '\t',
                               m, '\t',
                               len( row['seen']), '\t',
                               '%.2f' % (mean(row['num_concepts'])) if len(row['num_concepts'])>0 else 'null' , '\t',
                               '%.2f' % (mean(row['phi_concepts'])) if len(row['phi_concepts'])>0 else 'null' , '\t',
                               '%.2f' % (mean(row['phi_network'])) if len(row['phi_network'])>0 else 'null')

        if condition_name in ['1A', '1B']:
            for i, m in enumerate(motifs_mapping):
                row = seen[to_2d_tuple(normalize_cm(m))]
                print(*printer(i, m, row))
                
        else:
            for i, m in enumerate(motifs_plus_self_mapping):
                row = seen[to_2d_tuple(normalize_cm(m))]
                print(*printer(i, m, row))
            
            # for i, keyval in enumerate(seen.items()):
            #     key, row = keyval
            #     print(*printer(i, key, row))
    





    







#########################################################################################
# Utility Functions
#########################################################################################

def to_2d_tuple(arr):
    ''' takes a (mutable) 2D array and converts it to an (immutable) tuple '''
    return tuple(tuple(row) for row in arr)

def to_2d_list(arr):
    ''' takes a (immutable) 2D tuple and converts it to an (mutable) list '''
    return list(list(row) for row in arr)
    
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
    # return true if all nodes of incoming sum==x also have a mech in mechs
    
    in_count = count_incoming(cm)
    for count, node in zip(in_count, nodes):
        if (count == sm):
            if (node in mechs):
                continue
            return False
    return True
    

def gen_unique_cms(num_nodes):
    cms = product(*[[0,1]]*(num_nodes**2)) # Warning: mind-numbingly inefficient (brute force), don't use a lot of nodes (like, more than 4 sounds scary but I haven't tried it)

    unique_cms = set() #there will be 104 unique cms at 3 nodes with recurrent loops

    # look at all cms, add uniques to unique_cms
    for cm in cms:
        np_cm = np.array(cm)
        np_cm.shape = (num_nodes, num_nodes)
        unique_cms.add(tuple([tuple(x) for x in normalize_cm(np_cm)]))

    return unique_cms



if __name__ == "__main__":
    main()
