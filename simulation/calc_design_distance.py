import sys
import os

sys.path.append('./src')

from sage.all import *
from TamoBergTwoSets import *
import csv

q = 256  # Field size
n = 256  # Code dimension
r = [14, 3]  # Locality of the code
local_minimum_distance = [3, 3]  # correct one error
sub_group_type = ["add", "mult"]

# GF
F = GF(q, repr='int')

tmp = TamoBergCodeTwoSets(F, n, 1, r, local_minimum_distance, sub_group_type)
k = tmp.max_dimension()  # Code Dimension set to max possible.
C = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type)

sub_group, sub_group_type = C.sub_group()

## Simulation name and directory
sim_name = "GF_" + str(q) + "_" + sub_group_type[0] + "_" + str(len(sub_group[0])) + "_" + sub_group_type[1] + "_" + str(len(sub_group[1]))
res_dir = './results/' + sim_name

# Create Directory of not exists
isExist = os.path.exists(res_dir)
if not isExist:
    os.makedirs(res_dir)

res_path = res_dir + '/' + sim_name + '_design_distance.csv'
res_file_handle = open(res_path, 'w', newline='')

writer = csv.writer(res_file_handle)

print("## Setup Done ##")
print("## Start Calculation ", sim_name, "##")

max_k = tmp.max_dimension()
print("Max K:",max_k)
for k in range(1, max_k+1):
    tmp = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance,sub_group_type)
    d_lb = tmp.design_distance()
    print(k, ":", d_lb)
    writer.writerow([k, d_lb])

print("## Calculation Done ##")
res_file_handle.close()