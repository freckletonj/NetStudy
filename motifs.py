'''

'''

import numpy as np
from itertools import combinations, permutations, product
from normalize_cm import normalize_cm


  
# GENERATE TEST CMs #############################################################################
# number of unique_cms:
# 3 nodes = 104 : calculates fast
# 4 nodes = 3044 :calculates in <1 min

num_nodes = 4

cms = product(*[[0,1]]*(num_nodes**2)) # Warning: mind-numbingly inefficient, don't use a lot of nodes (like, more than 4 sounds scary but I haven't tried it)

unique_cms = set() #there will be 104 unique cms at 3 nodes with recurrent loops

for cm in cms:
    np_cm = np.array(cm)
    np_cm.shape = (num_nodes, num_nodes)
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
for cm in cms1:
    print(np.array(cm))



# 1.b) same as (1a) with the possibility of XORs for those nodes that
#       actually have 2 inputs





########################################################################################
# Condition (2)  -   All the nodes have self loops

# 2.a) only linear threshold mechanisms: >= 1, >= 2, >= 3


# 2.b) same as (2a) but with the possibility of XORs

     
