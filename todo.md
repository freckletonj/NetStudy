In this little project we wanted to check the
- average number of concepts
- average sum of the concepts' small phi values
- average big phi (probably always 0 for most of the motifs, though)
of all the 13 possible 3 node motifs.

for the following 4 conditions:

1. No self loops (that is each node has at most 2 inputs and 2 outputs):
     a) only linear threshold mechanisms: two possibilities Threshold >=
1 (Copy/OR), Threshold >= 2 (AND)
         maybe only those nodes that actually have 2 inputs should be
ANDs (not sure about this one)
     b) same as a) with the possibility of XORs for those nodes that
actually have 2 inputs

2. All the nodes have self loops
     a) only linear threshold mechanisms: >= 1, >= 2, >= 3
     b) same as a) but with the possibility of XORs

Could be plotted as bar graphs labeled by the motif number as in the
Sporns & Kötter 2004 paper Fig. 1.


So this is very similar to what you did for all possible 4 node
networks, just that the network structure is constrained by each motif
and then the average should be for each motif. So there are much less
possibilities of different mechanisms for these motifs than for the 4
node networks.

So for example motif 4 ( A <-> B <- C) could under condition 1a) have B
either an AND or an OR, but A must be a Copy and C is actually a "Null"
node, so it doesn't have inputs. We'll assume that Null nodes always
switch off.

Let me know if you have questions and if you would like to start I'd
suggest you try condition 1a) and then talk to me about the results
before doing the rest.
Best,

Larissa
