import sys
import os

sys.path.append('./src')

import csv
from networkx.algorithms import *
import matplotlib.pyplot as plt
import networkx as nx
from sage.coding.channel import StaticErrorRateChannel
from sage.coding.channel import random_error_vector
from TamoBargTwoSets import TamoBargCodeTwoSets
from sage.all import *
from Expander import *
from AdditiveSubgroup import *
from random import sample

# Simulation Control
####################

# "GlobalDecoder"
# "GlobalErasureErrorDecoder"
# "IterativeDecoder"
# "IterativeErasureErrorDecoder"
# "TwoStepsDecoder"
# "TwoStepsErasureErrorDecoder"
decoder_str = "IterativeDecoder"

# encode or use zero cw
use_enc = False

# excat distnace
exact_distance = None

sim_itr = 10  # Statistical Accuracy
print_freq_factor = 5  # Print frequency

start_err = 1

# Number of error to simulate
if exact_distance is None:
    sim_num_of_err = 16
else:
    sim_num_of_err = exact_distance + 1

q = 64  # Field size
n = 64  # Code Length
k = 6   # Code Dimension
r = [2, 2]  # Locality of the code
local_minimum_distance = [3, 3]  # correct one error
sub_group_type = ["add", "add"]
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

C = TamoBargCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type, shift_add=False, subgroup=add_subgroup)

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
res_dir = './results/' + sim_name + decoder_str

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

writer.writerow(["Num_of_error_symbols", \
    "probability_of_error_increase", "num_of_error_increase", \
    "probability_of_error_no_change","num_of_error_no_change", \
    "probability_of_error_decrease","num_of_error_decrease", \
    "probability_of_failure", "num_of_failure", \
    "avg_remin_err_weight", "min_remin_err_weight", "max_remin_err_weight"])

## Error Injection
#error_positions = get_error_position_subgraph(G,3)
#error_positions = [0, 4, 8, 12, 18, 22, 26, 30, 32, 36, 40, 44, 50, 54, 58, 62]
error_positions = list(range(0,n))
print("error_positions:", error_positions)

max_num_of_err = min(n-k, sim_num_of_err, len(error_positions))
print_freq = int(sim_itr / print_freq_factor)

for n_err in range(start_err, max_num_of_err + 1):
    print("Error Weight: ", n_err)
    twostep_fail_itr = 0
    iter_error_incr_itr = 0
    iter_error_decr_itr = 0
    iter_error_no_change = 0
    overall_remin_err_weight = 0
    max_remin_err_weight = 0
    min_remin_err_weight = n

    for i in range(0, sim_itr):
        error_positions_sample = sample(error_positions, n_err)
        err_vec = random_error_vector(n, F, error_positions_sample)
        r = c + err_vec
        if shorten:
            r[0] = c[0]

        e = r - c
        correct_c, num_of_itr = Dec.decode_to_code(r)
        remin_err_weight = (c - correct_c).hamming_weight()

        if remin_err_weight > max_remin_err_weight:
            max_remin_err_weight = remin_err_weight

        if remin_err_weight < min_remin_err_weight:
            min_remin_err_weight = remin_err_weight

        overall_remin_err_weight += remin_err_weight

        if remin_err_weight > n_err:
            iter_error_incr_itr += 1
        elif remin_err_weight < n_err:
            iter_error_decr_itr += 1
        else:
            iter_error_no_change += 1

        if n_err <=  ((C.design_distance() -1) // 2):
            if remin_err_weight > ((C.design_distance() -1) // 2):
                twostep_fail_itr += 1

        if (i + 1) % print_freq == 0:
            print("Iteration ", i + 1, "/", sim_itr, " Done.")

    writer.writerow([n_err, \
        100 * (iter_error_incr_itr / sim_itr) , iter_error_incr_itr, \
        100 * (iter_error_no_change / sim_itr) , iter_error_no_change, \
        100 * (iter_error_decr_itr / sim_itr) , iter_error_decr_itr, \
        100 * (twostep_fail_itr / sim_itr) , twostep_fail_itr, \
        (overall_remin_err_weight / sim_itr), min_remin_err_weight, max_remin_err_weight])

print("## Simulation Done ##")
res_file_handle.close()