'''

For Studying : http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0020369
And comparing their results against the theory of Integrated Information


'''

import numpy as np
from itertools import combinations, permutations, product
from normalize_cm import normalize_cm
from normalize_tpm import normalize_tpm
from generate_tpm import generate_tpm, nodes_to_short, AND, OR, XOR


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

    all_activations =  [p for p in product(*[[0,1]*NUM_NODES])] # all possible activation states


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

    has_self_loops = lambda cm: bool(len([1 for i in range(len(cm[0])) if cm[i][i] == 1]))
    has_unconnected_node = lambda cm: bool(len(
        [1 for i in range(len(cm[0])) if
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

    # Testing: Short Circuit these lists:
    #all_activations = [[0,1,1]]
    MECHANISMS = [OR, AND]

    duplicate_tpms = 0
    
    for cm in cms1:
        #print("CONNECTIVITY MATRIX: ", cm)

        seen_tpms_by_condition[cm] = set()

        for current_state in all_activations:
            #print("CURRENT STATE: ", current_state)

            for nodes in product(MECHANISMS, repeat=len(cm[0])):
                #print("NODES: ", nodes_to_short(nodes))

                tpm = normalize_tpm(generate_tpm(nodes, cm))

                # continue next loop if this tpm has already been seen on this condition
                if to_2d_tuple(tpm) in seen_tpms_by_condition[to_2d_tuple(cm)]:                    
                    duplicate_tpms += 1
                    continue

                seen_tpms_by_condition[cm].add(to_2d_tuple(tpm))

                # network = pyphi.Network(tpm, current_state, past_state, connectivity_matrix=cm)

                # subsystem = pyphi.Subsystem(range(network.size), network)

                # main_complex = pyphi.compute.main_complex(network)

    print("\n"*3)
    for key, value in seen_tpms_by_condition.items():
        print(key,"--", len(value))
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
    return tuple(tuple(row) for row in arr)




if __name__ == "__main__":
    main()
