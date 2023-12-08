import sys
import os

sys.path.append('./src')

import csv
from networkx.algorithms import *
import matplotlib.pyplot as plt
import networkx as nx
from sage.coding.channel import StaticErrorRateChannel
from TamoBergTwoSets import TamoBergCodeTwoSets
from sage.all import *
from AdditiveSubgroup import *

# Simulation Control
####################

# "GlobalDecoder"
# "GlobalErasureErrorDecoder"
# "IterativeDecoder"
# "IterativeErasureErrorDecoder"
# "TwoStepsDecoder"
# "TwoStepsErasureErrorDecoder"
decoder_str = "TwoStepsDecoder"

# encode or use zero cw
use_enc = False

# excat distnace
exact_distance = None

sim_itr = 5000  # Statistical Accuracy
print_freq_factor = 5  # Print frequency

# Number of error to simulate
if exact_distance is None:
    sim_num_of_err = 125
else:
    sim_num_of_err = exact_distance + 1

q = 64  # Field size
n = 64  # Code Length
k = 6   # Code Dimension
r = [2, 4]  # Locality of the code
local_minimum_distance = [3, 4]  # correct one error
sub_group_type = ["add", "mult"]
max_num_of_itr = 10

# Shorten code in case of different sub-group types.
if sub_group_type[0] != sub_group_type[1]:
    shorten = True
else:
    shorten = False

# GF
F = GF(q, repr='int')

# Specify Additive subgroup
if sub_group_type[0] == "add":
    sub_group_size = r[0] + local_minimum_distance[0] -1
    additive_subgroups = find_additive_subgroups(F, sub_group_size)

add_subgroup = [None, None]

C = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type, shift_add=False, subgroup=add_subgroup)

# Message Space
M = VectorSpace(F, k)
# Choose random message
message = M.random_element()

# Code space
V = VectorSpace(F, n)

Dec = C.decoder(decoder_str, max_num_of_itr)
Enc = C.encoder("VectorEncoder")

# Use Encoder or Zero Codeword
if use_enc:
    c = Enc.encode(message)
else:
    c = vector(F, [F.zero()] * n)  # All Zero CW

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

# Simulation name and directory
sim_name = "GF_" + str(q) + "_" + sub_group_type[0] + "_" + str(len(sub_group[0])) + "_" + sub_group_type[1] + "_" + str(len(sub_group[1]))
res_dir = './results/' + sim_name + "_" + decoder_str

# Create Directory of not exists
isExist = os.path.exists(res_dir)
if not isExist:
    os.makedirs(res_dir)

res_path = res_dir + '/' + sim_name + '.csv'
graph_path = res_dir + '/' + sim_name + '_graph.png'
res_file_handle = open(res_path, 'w', newline='')

writer = csv.writer(res_file_handle)

print("## Setup Done ##")
print("## Start Simulation ", sim_name, "##")
print(C)
print("GF: ", q)
print("Code dim n: ", n)
print("Message dim k: ", k)
print("k_max: ", C.max_dimension())
print("max rate: ", C.max_rate())
print("Locality r1: ", r[0])
print("Locality r2: ", r[1])
print("Local Minimum distance d1: ", local_minimum_distance[0])
print("Local Minimum distance d2: ", local_minimum_distance[1])
print("Design Distance: ", C.design_distance())
print("First Subgroup: ", sub_group[0])
print("First Subgroup Type: ", sub_group_type[0])
print("First Subgroup Size: ", len(sub_group[0]))
print("Second Subgroup: ", sub_group[1])
print("Second Subgroup Type: ", sub_group_type[1])
print("Second Subgroup Size: ", len(sub_group[1]))
print("Partitions Graph is Connected:", nx.is_connected(G))
print("Number of disjoint subgraphs:", nx.number_connected_components(G))
print("Evaluation Points: ", evalpts)
print("Codeword : ", c)

# Save to CSV
writer.writerow(["GF: ", q])
writer.writerow(["Code dim n: ", n])
writer.writerow(["Message dim k: ", k])
writer.writerow(["k_max: ", C.max_dimension()])
writer.writerow(["max rate: ", C.max_rate()])
writer.writerow(["Locality r1: ", r[0]])
writer.writerow(["Locality r2: ", r[1]])
writer.writerow(["Local Minimum distance d1: ", local_minimum_distance[0]])
writer.writerow(["Local Minimum distance d2: ", local_minimum_distance[1]])
writer.writerow(["Design Distance: ", C.design_distance()])
writer.writerow(["First Subgroup: ", sub_group[0]])
writer.writerow(["First Subgroup Type: ", sub_group_type[0]])
writer.writerow(["First Subgroup Size: ", len(sub_group[0])])
writer.writerow(["Second Subgroup: ", sub_group[1]])
writer.writerow(["Second Subgroup Type: ", sub_group_type[1]])
writer.writerow(["Second Subgroup Size: ", len(sub_group[1])])
writer.writerow(["Evaluation Points: ", evalpts])
writer.writerow(["Partitions Graph is Connected:", nx.is_connected(G)])
writer.writerow(["Number of disjoint subgraphs:",nx.number_connected_components(G)])
writer.writerow(["Partitions Graph Edge Expansion:",cuts.edge_expansion(G, G_R)])
writer.writerow(["Codeword : ", c])

# Save Graph Picture
plt.savefig(graph_path)

max_num_of_err = min(n-k, sim_num_of_err)
print_freq = int(sim_itr / print_freq_factor)

writer.writerow(["Num_of_error_symbols", "average_num_of_iteration", "probability_of_success"])

for n_err in range(1, max_num_of_err):
    print("Error Weight: ", n_err)
    Chan = StaticErrorRateChannel(V, n_err)
    success_itr = 0
    overall_num_itr = 0

    for i in range(0, sim_itr):
        r = Chan.transmit(c)
        if shorten:
            r[0] = c[0]

        e = r - c
        if(decoder_str == "TwoStepsErasureErrorDecoder"):
            correct_c, num_of_itr, _ = Dec.decode_to_code(r)
        else:
            correct_c, num_of_itr = Dec.decode_to_code(r)

        if(decoder_str == "IterativeErasureErrorDecoder"):
            correct_c = correct_c[0]

        if (correct_c == c):
            success = True
            success_itr += 1
            overall_num_itr += num_of_itr
        else:
            success = False

        if (i + 1) % print_freq == 0:
            print("Iteration ", i + 1, "/", sim_itr, " Done.")

    if success_itr != 0:
        avg_num_of_itr = (overall_num_itr / success_itr)
    else:
        avg_num_of_itr = max_num_of_itr

    writer.writerow([n_err, avg_num_of_itr, 100 * (success_itr / sim_itr)])

print("## Simulation Done ##")
res_file_handle.close()