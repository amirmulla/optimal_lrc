import sys
from numpy import double
sys.path.append('./src')

from sage.all import *
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
from networkx.algorithms import *
from TamoBergTwoSets import TamoBergCodeTwoSets
from Expander import *
from AdditiveSubgroup import *
from itertools import combinations
import math

def findsubsets(s, n):
    return list(map(set, combinations(s, n)))

#set_random_seed(100)


q = 64  # Field size
n = 64  # Code dimension
k = 8  # Information/message dimension
r = [6, 2]  # Locality of the code
local_minimum_distance = [3, 2]  # correct one error
sub_group_type = ["add", "mult"]
n_err = 2
max_num_of_itr = 10

print("Factor q:", factor(q))
print("Factor q-1:", factor(q-1))

# GF
F = GF(q, name="a", repr='log')

# Code space
V = VectorSpace(F, n)

# Message Space
M = VectorSpace(F, k)
# Choose random message
message = M.random_element()

#å Specify Additive subgroup
if sub_group_type[0] == "add":
    sub_group_size = r[0] + local_minimum_distance[0]-1
    additive_subgroups = find_additive_subgroups(F, sub_group_size)

#add_subgroup = [additive_subgroups[0], None]
add_subgroup = [None, None]

#C = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type, shift_add=True, subgroup=add_subgroup)
C = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type)

print("GF: ", q)
print("Code dim n: ", n)
print("Message dim k: ", k)
print("Locality r1: ", r[0])
print("Locality r2: ", r[1])

sub_group, sub_group_type = C.sub_group()

G, G_L, G_R = C.bipartite_graph()
m, n = len(G_L), len(G_R)
pos = dict()
pos.update((i, (i - m / 2, 1)) for i in range(m))
pos.update((i, (i - m - n / 2, 0)) for i in range(m, m + n))
nx.draw(G, with_labels=True, pos=pos, node_size=300, width=0.4)
plt.savefig('graph.png')

print("Partitions Graph is Connected:", nx.is_connected(G))
print("Number of disjoint subgraphs:", nx.number_connected_components(G))

print("First Subgroup: ", sub_group[0])
print("First Subgroup Type: ", sub_group_type[0])
print("First Subgroup Size: ", len(sub_group[0]))
print("Second Subgroup: ", sub_group[1])
print("Second Subgroup Type: ", sub_group_type[1])
print("Second Subgroup Size: ", len(sub_group[1]))

print(G)
print(G_L)
print(len(G_L))
print(G_R)
print(len(G_R))
print("##########################")

if len(G_L) > len(G_R):
    S = G_L
    d = len(sub_group[0])
    if sub_group_type[0] != sub_group[1]:
        S.remove(0)
else:
    S = G_R
    d = len(sub_group[1])

print("S:",S)
print("##########################")


#print(nx.normalized_laplacian_spectrum(G))


l_adj = nx.adjacency_spectrum(G)
l_new = []

for l in l_adj:
    l_new.append(round(float(l),2))

l_new.sort()
print("Eigen values: ",l_new)
print("Spectral Gap: ",l_new[-1] - l_new[-2])
print("algebric connectivity: ",nx.algebraic_connectivity(G=G, method="tracemin_lu"))

print("Adj matrix:")
print(nx.to_pandas_adjacency(G))


eps = float(1/d)

delta = d # max delta value
for i in range(0, int(math.floor(eps*len(S)))):
    combinations_list = findsubsets(S, i+1)
    for comb in combinations_list:
        tmp = nx.node_expansion(G, comb)
        if tmp < delta:
            delta = tmp

    print(i+1, "/",int(eps*len(S)),":", delta)
