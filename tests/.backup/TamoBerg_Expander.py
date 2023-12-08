from cProfile import label
import sys
sys.path.append('./src')

from sage.all import *
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
from networkx.algorithms import *
from networkx.algorithms.connectivity import EdgeComponentAuxGraph
from sage.coding.channel import StaticErrorRateChannel
from TamoBergTwoSets import TamoBergCodeTwoSets
from Expander import *
from AdditiveSubgroup import *


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
#message = zero_vector(F,k)

# Specify Additive subgroup
if sub_group_type[0] == "add":
    sub_group_size = r[0] + local_minimum_distance[0]-1
    additive_subgroups = find_additive_subgroups(F, sub_group_size)

add_subgroup = [None, None]

C = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type, shift_add=False, subgroup=add_subgroup)

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


def nudge(pos, x_shift, y_shift):
    return {n: (x + x_shift, y + y_shift) for n, (x, y) in pos.items()}


G, G_L, G_R = C.bipartite_graph()
m, n = len(G_L), len(G_R)
pos = dict()
pos.update((i, (1,i - m/2)) for i in range(m))
pos.update((i, (2,i - m - n/2)) for i in range(m, m + n))


fig, ax = plt.subplots()

labels = nx.get_node_attributes(G, 'coset')
#pos_nodes = nudge(pos, 0, 0.2)                                # shift the layout
nx.draw_networkx(G, pos=pos, with_labels=False, ax=ax)         # default nodes and edges

l_pos_nodes = dict()
r_pos_nodes = dict()
l_pos_nodes.update((i, (1, i - m/2)) for i in range(m))
l_pos_nodes = nudge(l_pos_nodes, x_shift=-0.25, y_shift=0)
r_pos_nodes.update((i, (2, i - m - n/2)) for i in range(m, m + n))
r_pos_nodes = nudge(r_pos_nodes, x_shift=0.15, y_shift=0)

pos_nodes = dict()
pos_nodes.update(l_pos_nodes)
pos_nodes.update(r_pos_nodes)
nx.draw_networkx_labels(G, pos=pos_nodes, ax=ax, labels=labels)  # nudged labels
#ax.set_xlim(tuple(i*1.1 for i in ax.get_xlim()))
print(ax.get_xlim())
ax.set_xlim(0.5, 2.5)
plt.show()
#plt.savefig('graph.png')


exit()
print("Partitions Graph is Connected:", nx.is_connected(G))
print("Number of disjoint subgraphs:", nx.number_connected_components(G))

for i in range(0,nx.number_connected_components(G)):
    print(get_error_position_subgraph(G,i))

#G_L = {n for n, d in G.nodes(data=True) if d["bipartite"] == 0}
#nx.draw(G.subgraph(x), with_labels=True, pos=pos, node_size=300, width=0.4)
#plt.show()
