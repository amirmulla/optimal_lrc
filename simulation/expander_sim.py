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

sim_itr = 10000  # Statistical Accuracy
print_freq_factor = 5  # Print frequency

start_err = 1

# Number of error to simulate
if exact_distance is None:
    sim_num_of_err = 64
else:
    sim_num_of_err = exact_distance + 1

q = 256  # Field size
n = 256  # Code dimension
k = 6  # Information/message dimension
r = [4, 4]  # Locality of the code
local_minimum_distance = [5, 5]  # correct one error
sub_group_type = ["add", "add"]
max_num_of_itr = 10
biased = False
t = 64 # error boundary [0,t]

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
if biased == True:
    res_dir = './results/' + sim_name + "_expander_biased_" + decoder_str
else:
    res_dir = './results/' + sim_name + "_expander_uniform_" + decoder_str


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
if biased == True:
    error_positions = list(range(0,t))
else:
    error_positions = list(range(0,n))

# 81 (4x5)
#error_positions = [16, 21, 26, 28, 29, 31, 38, 44, 45, 49, 52, 53, 54, 59, 60, 62, 64, 66, 69, 77]
#error_positions = [2, 4, 5, 6, 10, 11, 18, 23, 30, 34, 41, 43, 47, 48, 57, 61, 68, 70, 74, 75]
#error_positions = [0, 1, 3, 7, 13, 15, 22, 24, 27, 33, 36, 42, 46, 50, 55, 58, 63, 65, 73, 76]
#error_positions = [8, 9, 12, 14, 17, 19, 20, 25, 32, 35, 37, 39, 40, 51, 56, 67, 71, 72, 78, 79]

# 256 (8x8)
#error_positions = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63]
#error_positions = [192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254, 255]
#error_positions = [128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191]
#error_positions = [64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127]

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