import sys
import os

sys.path.append('./src')

from sage.all import *
from TamoBergTwoSets import *
from sage.coding.channel import StaticErrorRateChannel

import csv

print_log = True

q = 256  # Field size
n = 256  # Code dimension
k = 230  # Information/message dimension
r = [14, 14]  # Locality of the code
local_minimum_distance = [3, 3]  # correct one error
sub_group_type = ["add", "add"]
max_num_of_itr = 10

# Puncture code
puncture = False

# GF
F = GF(q, repr='int')

# Code space
V = VectorSpace(F, n)

C = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type)
Dec = C.decoder("IterativeDecoder", max_num_of_itr)
# All zero code word
c = vector(F, [F.zero()] * n)

sub_group, sub_group_type = C.sub_group()
partitions, _ = C.partition()
evalpts = C.evaluation_points()

sim_name = "GF_" + str(q) + "_" + sub_group_type[0] + "_" + str(len(sub_group[0])) + "_" + sub_group_type[1] + "_" + str(len(sub_group[1]))
res_dir = './results'

# Create Directory of not exists
isExist = os.path.exists(res_dir)
if not isExist:
    os.makedirs(res_dir)

res_path = res_dir + '/' + sim_name + '.csv'
log_path = res_dir + '/' + sim_name + '_log.txt'
res_file_handle = open(res_path, 'w', newline='')
log_file_handle = open(log_path, 'w', newline='')

writer = csv.writer(res_file_handle)

print("## Setup Done ##")
print("## Start Simulation ", sim_name, "##")
print(C)

print("GF: ", q)
print("Code dim n: ", n)
print("Message dim k: ", k)
print("Locality r1: ", r[0])
print("Locality r2: ", r[1])
print("Local Minimum distance d1: ", local_minimum_distance[0])
print("Local Minimum distance d2: ",local_minimum_distance[1])
print("First Subgroup: ", sub_group[0])
print("First Subgroup Type: ", sub_group_type[0])
print("First Subgroup Size: ", len(sub_group[0]))
print("Second Subgroup: ", sub_group[1])
print("Second Subgroup Type: ", sub_group_type[1])
print("Second Subgroup Size: ", len(sub_group[1]))
print("Evaluation Points: ", evalpts)

# Save to CSV
writer.writerow(["GF: ", q])
writer.writerow(["Code dim n: ", n])
writer.writerow(["Message dim k: ", k])
writer.writerow(["Locality r1: ", r[0]])
writer.writerow(["Locality r2: ", r[1]])
writer.writerow(["Local Minimum distance d1: ", local_minimum_distance[0]])
writer.writerow(["Local Minimum distance d2: ", local_minimum_distance[1]])
writer.writerow(["First Subgroup: ", sub_group[0]])
writer.writerow(["First Subgroup Type: ", sub_group_type[0]])
writer.writerow(["First Subgroup Size: ", len(sub_group[0])])
writer.writerow(["Second Subgroup: ", sub_group[1]])
writer.writerow(["Second Subgroup Type: ", sub_group_type[1]])
writer.writerow(["Second Subgroup Size: ", len(sub_group[1])])
writer.writerow(["Evaluation Points: ", evalpts])

if print_log:
    print(C, file=log_file_handle)
    print("GF: ", q, file=log_file_handle)
    print("Code dim n: ", n, file=log_file_handle)
    print("Message dim k: ", k, file=log_file_handle)
    print("Locality r1: ", r[0], file=log_file_handle)
    print("Locality r2: ", r[1], file=log_file_handle)
    print("Local Minimum distance d1: ", local_minimum_distance[0], file=log_file_handle)
    print("Local Minimum distance d2: ",local_minimum_distance[1], file=log_file_handle)
    print("First Subgroup: ", sub_group[0], file=log_file_handle)
    print("First Subgroup Type: ", sub_group_type[0], file=log_file_handle)
    print("First Subgroup Size: ", len(sub_group[0]), file=log_file_handle)
    print("Second Subgroup: ", sub_group[1], file=log_file_handle)
    print("Second Subgroup Type: ", sub_group_type[1], file=log_file_handle)
    print("Second Subgroup Size: ", len(sub_group[1]), file=log_file_handle)
    print("Evaluation Points: ", evalpts, file=log_file_handle)

max_num_of_err = 32
sim_itr = 5000

print_freq = int(sim_itr/5)

writer.writerow(["Num_of_error_symbols", "average_num_of_iteration","probability_of_success"])

for n_err in range(1, max_num_of_err):
    print("Error Weight: ", n_err)
    Chan = StaticErrorRateChannel(V, n_err)
    success_itr = 0
    overall_num_itr = 0

    for i in range(0, sim_itr):
        r = Chan.transmit(c)
        if puncture:
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
            overall_num_itr += num_of_itr
        else:
            success = False

        if (i+1) % print_freq == 0:
            print("Iteration ", i+1, "/", sim_itr," Done.")

    if success_itr != 0:
        avg_num_of_itr = (overall_num_itr / success_itr)
    else:
        avg_num_of_itr = max_num_of_itr

    writer.writerow([n_err, avg_num_of_itr, 100 * (success_itr / sim_itr)])

    if print_log:
        print("Num of error symbols: ", n_err, file=log_file_handle)
        print("Average Num of Iteration: ", avg_num_of_itr,file=log_file_handle)
        print("Probability of Success: ",100 * (success_itr / sim_itr),file=log_file_handle)
        print("##########################", file=log_file_handle)

print("## Simulation Done ##")
log_file_handle.close()
res_file_handle.close()