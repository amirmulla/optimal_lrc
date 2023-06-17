import sys
import os
sys.path.append('./src')

from TamoBergTwoSets import TamoBergCodeTwoSets

import csv
from sage.all import *
from AdditiveSubgroup import *
import matplotlib.pyplot as plt
from networkx.algorithms import *
import networkx as nx

q = 64  # Field size
n = 64  # Code dimension
k = 6  # Information/message dimension
r = [2, 5]  # Locality of the code
local_minimum_distance = [3, 3]  # correct two errors
sub_group_type = ["add", "mult"]

# GF
F = GF(q, repr='int')

# Specify Additive subgroup
if sub_group_type[0] == "add":
    sub_group_size = r[0] + local_minimum_distance[0] - 1
    additive_subgroups = find_additive_subgroups(F, sub_group_size)

# GF(64) 4x7 --- idx 10
add_idx = None

if add_idx is None:   
    add_subgroup = [None, None]
else:
    add_subgroup = [additive_subgroups[add_idx], None]

C = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type, shift_add=False, subgroup=add_subgroup)

sub_group, sub_group_type = C.sub_group()
partitions, _ = C.partition()
evalpts = C.evaluation_points()

# Bipartite Graph
G, G_L, G_R = C.bipartite_graph()
left, right = len(G_L), len(G_R)
pos = dict()
pos.update((i, (i - left / 2, 1)) for i in range(left))
pos.update((i, (i - left - right / 2, 0)) for i in range(left, left + right))
nx.draw(G, with_labels=True, pos=pos, node_size=300, width=0.4)

print(C)
print("k_max: ", C.max_dimension())
print("Design Distance: ", C.design_distance())
print("Local Minimum distance d1: ", local_minimum_distance[0])
print("Local Minimum distance d2: ", local_minimum_distance[1])
print("First Subgroup: ", sub_group[0])
print("First Subgroup Type: ", sub_group_type[0])
print("First Subgroup Size: ", len(sub_group[0]))
print("Second Subgroup: ", sub_group[1])
print("Second Subgroup Type: ", sub_group_type[1])
print("Second Subgroup Size: ", len(sub_group[1]))
print("Partitions Graph is Connected:", nx.is_connected(G))
print("Number of disjoint subgraphs:", nx.number_connected_components(G))
print("Partitions Graph Edge Expansion:", cuts.edge_expansion(G, G_R))

if add_idx is None:
    sim_name = "GF_" + str(q) + "_" + sub_group_type[0] + "_" + str(len(sub_group[0])) + "_" + sub_group_type[1] + "_" + str(len(sub_group[1]))
else:    
    sim_name = "GF_" + str(q) + "_" + sub_group_type[0] + "_" + str(len(sub_group[0])) + "_" + sub_group_type[1] + "_" + str(len(sub_group[1])) + "_idx_" + str(add_idx)

res_dir = './results/' + sim_name + "_weight_dist"

print(sim_name)

# Create Directory of not exists
isExist = os.path.exists(res_dir)
if not isExist:
    os.makedirs(res_dir)

res_path = res_dir + '/' + sim_name + '_weight_dist.csv'
res_file_handle = open(res_path, 'w', newline='')
writer = csv.writer(res_file_handle)

writer.writerow(["k_max: ", C.max_dimension()])
writer.writerow(["Design Distance: ", C.design_distance()])
writer.writerow(["Partitions Graph is Connected:", nx.is_connected(G)])
writer.writerow(["Number of disjoint subgraphs:", nx.number_connected_components(G)])
writer.writerow(["Number of disjoint subgraphs:", nx.number_connected_components(G)])
writer.writerow(["Partitions Graph Edge Expansion:", cuts.edge_expansion(G, G_R)])

print("## Setup Done ##")
print("## Start Simulation ", sim_name, "##")
w = C.weight_distribution()
writer.writerow(["Weight Distribution: ", w])

print(w)

for i in range(1,len(w)):
    if w[i] != 0:
        print("Code Distance: ", i)
        writer.writerow(["Code Distance: ", i])
        break

print("## Simulation Done ##")
res_file_handle.close()
