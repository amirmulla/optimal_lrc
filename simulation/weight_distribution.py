import sys
import os
sys.path.append('./src')

from TamoBergTwoSets import TamoBergCodeTwoSets

import csv
from sage.all import *

q = 25  # Field size
n = 24  # Code Length
k = 7  # Code Dimension
r = [2, 2]  # Locality of the code
local_minimum_distance = [3, 2]  # correct one error
sub_group_type = ["mult", "mult"]

# GF
F = GF(q, repr='int')

C = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type, shift_add=False)
sub_group, sub_group_type = C.sub_group()
partitions, _ = C.partition()
evalpts = C.evaluation_points()

print(C)
print("k_max: ", C.max_dimension())
print("Design Distance: ", C.design_distance())

sim_name = "GF_" + str(q) + "_" + sub_group_type[0] + "_" + str(len(sub_group[0])) + "_" + sub_group_type[1] + "_" + str(len(sub_group[1]))
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
