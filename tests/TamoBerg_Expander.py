import sys

sys.path.append('./src')

from sage.all import *
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
from networkx.algorithms import *
from networkx.algorithms.connectivity import EdgeComponentAuxGraph

from TamoBergTwoSets import *
from sage.coding.channel import StaticErrorRateChannel


#set_random_seed(100)

q = 64  # Field size
n = 64  # Code dimension
k = 16  # Information/message dimension
r = [2, 2]  # Locality of the code
local_minimum_distance = [3, 3]  # correct one error
sub_group_type = ["add", "add"]
n_err = 8
max_num_of_itr = 10

# GF
F = GF(q, name="a", repr='int')

# Code space
V = VectorSpace(F, n)

# Message Space
M = VectorSpace(F, k)
# Choose random message
message = M.random_element()
#message = zero_vector(F,k)

C = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type)

print(C)
sub_group, sub_group_type = C.sub_group()
partitions, _ = C.partition()
evalpts = C.evaluation_points()

print("partitions:",partitions)

G, G_L, G_R = C.bipartite_graph()
m, n = len(G_L), len(G_R)
pos = dict()
pos.update((i, (i - m / 2, 1)) for i in range(m))
pos.update((i, (i - m - n / 2, 0)) for i in range(m, m + n))
#nx.draw(G, with_labels=True, pos=pos, node_size=300, width=0.4)
#plt.savefig('graph.png')

print("Partitions Graph is Connected:", nx.is_connected(G))
print("Number of disjoint subgraphs:", nx.number_connected_components(G))

aux_graph = EdgeComponentAuxGraph.construct(G)
l_subgraphs = sorted(map(sorted, aux_graph.k_edge_components(k=1)))
x = set(l_subgraphs[1])
nx.draw(G.subgraph(x), with_labels=True, pos=pos, node_size=300, width=0.4)
plt.show()
