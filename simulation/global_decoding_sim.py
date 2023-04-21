import sys
import os

sys.path.append('./src')

from sage.all import *
from TamoBergTwoSets import *
from sage.coding.channel import StaticErrorRateChannel
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import bipartite
from networkx.algorithms import *
import csv


# Simulation Control
print_log = False  # Save log into file
use_erasure_decoder = False
use_enc = False
sim_itr = 5000  # Statistical Accuracy
print_freq_factor = 5  # Print frequency
sim_num_of_err = 64  # Number of error to simulate

q = 64  # Field size
n = 64  # Code Length
k = 10   # Code Dimension
r = [2, 4]  # Locality of the code
local_minimum_distance = [3, 5]  # correct one error
sub_group_type = ["add", "add"]
max_num_of_itr = 10

# Shorten code in case of different sub-group types.
if sub_group_type[0] != sub_group_type[1]:
    shorten = True
else:
    shorten = False

# GF
F = GF(q, repr='int')

# Code
C = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type)

# Message Space
M = VectorSpace(F, k)
# Choose random message
message = M.random_element()

# Code space
V = VectorSpace(F, n)

if use_erasure_decoder:
    Dec = C.decoder("IterativeErasureErrorDecoder", max_num_of_itr)
else:
    Dec = C.decoder("IterativeDecoder", max_num_of_itr)

Enc = C.encoder("VectorEncoder")

# Use Encoder or Zero Codeword
if use_enc:
    c = Enc.encode(message)
else:
    c = vector(F, [F.zero()] * n)  # All Zero CW

sub_group, sub_group_type = C.sub_group()
partitions, _ = C.partition()
evalpts = C.evaluation_points()

## Bipartite Graph
G, G_L, G_R = C.bipartite_graph()
left, right = len(G_L), len(G_R)
pos = dict()
pos.update((i, (i - left / 2, 1)) for i in range(left))
pos.update((i, (i - left - right / 2, 0)) for i in range(left, left + right))
nx.draw(G, with_labels=True, pos=pos, node_size=300, width=0.4)

## Simulation name and directory
sim_name = "GF_" + str(q) + "_" + sub_group_type[0] + "_" + str(
    len(sub_group[0])) + "_" + sub_group_type[1] + "_" + str(len(sub_group[1]))
if use_erasure_decoder:
    res_dir = './results/' + sim_name + "_erasure"
else:
    res_dir = './results/' + sim_name

# Create Directory of not exists
isExist = os.path.exists(res_dir)
if not isExist:
    os.makedirs(res_dir)

res_path = res_dir + '/' + sim_name + '.csv'
res_err_path = res_dir + '/' + sim_name + '_err.csv'

log_path = res_dir + '/' + sim_name + '_log.txt'
graph_path = res_dir + '/' + sim_name + '_graph.png'

res_file_handle = open(res_path, 'w', newline='')
res_err_file_handle = open(res_err_path, 'w', newline='')
log_file_handle = open(log_path, 'w', newline='')

writer = csv.writer(res_file_handle)
err_writer = csv.writer(res_err_file_handle)

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
print("Partitions Graph Edge Expansion:", cuts.edge_expansion(G, G_R))
print("Evaluation Points: ", evalpts)
print("Codeword : ", c)
print("Erasure Decoder: ", use_erasure_decoder)

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
writer.writerow(["Erasure Decoder: ", use_erasure_decoder])
# Save Graph Picture
plt.savefig(graph_path)

if print_log:
    print(C, file=log_file_handle)
    print("GF: ", q, file=log_file_handle)
    print("Code dim n: ", n, file=log_file_handle)
    print("Message dim k: ", k, file=log_file_handle)
    print("k_max: ", C.max_dimension(), file=log_file_handle)
    print("max rate: ", C.max_rate(), file=log_file_handle)
    print("Locality r1: ", r[0], file=log_file_handle)
    print("Locality r2: ", r[1], file=log_file_handle)
    print("Local Minimum distance d1: ",local_minimum_distance[0],file=log_file_handle)
    print("Local Minimum distance d2: ",local_minimum_distance[1],file=log_file_handle)
    print("Design Distance: ", C.design_distance(), file=log_file_handle)
    print("First Subgroup: ", sub_group[0], file=log_file_handle)
    print("First Subgroup Type: ", sub_group_type[0], file=log_file_handle)
    print("First Subgroup Size: ", len(sub_group[0]), file=log_file_handle)
    print("Second Subgroup: ", sub_group[1], file=log_file_handle)
    print("Second Subgroup Type: ", sub_group_type[1], file=log_file_handle)
    print("Second Subgroup Size: ", len(sub_group[1]), file=log_file_handle)
    print("Evaluation Points: ", evalpts, file=log_file_handle)
    print("Partitions Graph is Connected:",nx.is_connected(G),file=log_file_handle)
    print("Number of disjoint subgraphs:",nx.number_connected_components(G),file=log_file_handle)
    print("Partitions Graph Edge Expansion:",cuts.edge_expansion(G, G_R),file=log_file_handle)
    print("Codeword : ", c, file=log_file_handle)

max_num_of_err = sim_num_of_err
print_freq = int(sim_itr / print_freq_factor)

writer.writerow(["Num_of_error_symbols", "probability_of_success"])
err_writer.writerow(["Num_of_error_symbols", "remain_err_weight", "count", "dist"])

dd = C.design_distance()
dd_th = (dd - 1) // 2
for n_err in range(1,max_num_of_err):
    print("Error Weight: ", n_err)
    Chan = StaticErrorRateChannel(V, n_err)
    success_itr = 0
    overall_num_itr = 0
    remain_err_dist = [0] * (max_num_of_err + 4)

    for i in range(0, sim_itr):
        r = Chan.transmit(c)
        if shorten:
            r[0] = c[0]

        e = r - c
        correct_c, num_of_itr = Dec.decode_to_code(r)

        if print_log:
            print("Error              : ", e, file=log_file_handle)
            print("Corrected Codeword : ", correct_c, file=log_file_handle)
            print("Correction Successfull: ", correct_c == c, file=log_file_handle)

        if (correct_c == c):
            success = True
            success_itr += 1
            remin_err_weight = 0
        else:
            remin_err_weight = (c - correct_c).hamming_weight()
            if remin_err_weight <= dd_th: # still can be decoded by global decoding
                success = True
                success_itr += 1
            else:
                success = False

        remain_err_dist[remin_err_weight] += 1

        if (i + 1) % print_freq == 0:
            print("Iteration ", i + 1, "/", sim_itr, " Done.")

    writer.writerow([n_err, 100 * (success_itr / sim_itr)])

    for i in range(0,len(remain_err_dist)):
        if remain_err_dist[i] != 0:
            err_writer.writerow([n_err, i, remain_err_dist[i], remain_err_dist[i] / sim_itr])

    if print_log:
        print("Num of error symbols: ", n_err, file=log_file_handle)
        print("Probability of Success: ",100 * (success_itr / sim_itr),file=log_file_handle)
        print("##########################", file=log_file_handle)

print("## Simulation Done ##")

log_file_handle.close()
res_file_handle.close()
res_err_file_handle.close()