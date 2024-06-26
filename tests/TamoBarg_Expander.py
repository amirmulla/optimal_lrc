import sys
sys.path.append('./src')

from sage.all import *
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
from networkx.algorithms import *
from networkx.algorithms.connectivity import EdgeComponentAuxGraph
from sage.coding.channel import StaticErrorRateChannel
from TamoBargTwoSets import TamoBargCodeTwoSets
from Expander import *
from AdditiveSubgroup import *


#set_random_seed(100)

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

# Code space
V = VectorSpace(F, n)

# Message Space
M = VectorSpace(F, k)
# Choose random message
message = M.random_element()
#message = zero_vector(F,k)

# Specify Additive subgroup
if sub_group_type[0] == "add":
    sub_group_size = r[0] + local_minimum_distance[0]-1
    additive_subgroups = find_additive_subgroups(F, sub_group_size)

add_subgroup = [None, None]

C = TamoBargCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type, shift_add=False, subgroup=add_subgroup)

print("GF: ", q)
print("Code dim n: ", n)
print("Message dim k: ", k)
print("Max k: ", C.max_dimension())
print("Locality r1: ", r[0])
print("Locality r2: ", r[1])

sub_group, sub_group_type = C.sub_group()
partitions, _ = C.partition()
evalpts = C.evaluation_points()

print("Local Minimum distance d1: ", local_minimum_distance[0])
print("Local Minimum distance d2: ", local_minimum_distance[1])
print("Design Distance: ", C.design_distance())
print("First Subgroup: ", sub_group[0])
print("First Subgroup Type: ", sub_group_type[0])
print("First Subgroup Size: ", len(sub_group[0]))
print("Second Subgroup: ", sub_group[1])
print("Second Subgroup Type: ", sub_group_type[1])
print("Second Subgroup Size: ", len(sub_group[1]))

for p in partitions:
    i = 1
    for coset in p:
        print(i, ": ",coset)
        i += 1
    print("########")

G, G_L, G_R = C.bipartite_graph()
m, n = len(G_L), len(G_R)
pos = dict()
pos.update((i, (i - m / 2, 1)) for i in range(m))
pos.update((i, (i - m - n / 2, 0)) for i in range(m, m + n))
#nx.draw(G, with_labels=True, pos=pos, node_size=300, width=0.4)
#plt.savefig('graph.png')

print("Partitions Graph is Connected:", nx.is_connected(G))
print("Number of disjoint subgraphs:", nx.number_connected_components(G))

for i in range(0,nx.number_connected_components(G)):
    print(get_error_position_subgraph(G,i))

#G_L = {n for n, d in G.nodes(data=True) if d["bipartite"] == 0}
#nx.draw(G.subgraph(x), with_labels=True, pos=pos, node_size=300, width=0.4)
#plt.show()
