from sage.all import *
import sys
sys.path.append('./src')

from TamoBergTwoSets import TamoBergCodeTwoSets
from networkx.algorithms import *
import networkx as nx

q = 64  # Field size
n = 64  # Code dimension
k = 10  # Information/message dimension
r = [2, 2]  # Locality of the code
local_minimum_distance = [3, 3]  # correct one error
sub_group_type = ["add", "add"]
n_err = 2
max_num_of_itr = 10

print("Factor q:", factor(q))
print("Factor q-1:", factor(q-1))

# GF
F = GF(q, name="a", repr='int')

# Specify Additive subgroup
# if sub_group_type[0] == "add":
#    sub_group_size = r[0] + local_minimum_distance[0]-1
#    additive_subgroups = find_additive_subgroups(F, sub_group_size)
add_subgroup = [None, None]

C = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance,
                        sub_group_type, shift_add=False, subgroup=add_subgroup)

print("GF: ", q)
print("Code dim n: ", n)
print("Message dim k: ", k)
print("Locality r1: ", r[0])
print("Locality r2: ", r[1])

sub_group, sub_group_type = C.sub_group()

print("First Subgroup: ", sub_group[0])
print("First Subgroup Type: ", sub_group_type[0])
print("First Subgroup Size: ", len(sub_group[0]))
print("Second Subgroup: ", sub_group[1])
print("Second Subgroup Type: ", sub_group_type[1])
print("Second Subgroup Size: ", len(sub_group[1]))

# Spectral Analysis
G, G_L, G_R = C.bipartite_graph()
print("Partitions Graph is Connected:", nx.is_connected(G))
print("Number of disjoint subgraphs:", nx.number_connected_components(G))

l_adj = nx.adjacency_spectrum(G)
l_new = []

for l in l_adj:
    l_new.append(round(real(l), 2))

l_new.sort()
print("Eigen values: ", l_new)
print("Spectral Gap: ", l_new[-1] - l_new[-2])
print("algebric connectivity: ",
      nx.algebraic_connectivity(G=G, method="tracemin_lu"))

# print("Adj matrix:")
# print(nx.to_pandas_adjacency(G))
