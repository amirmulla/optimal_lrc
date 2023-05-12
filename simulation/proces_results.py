import csv
import os

q = 64  # Field size
n = 64  # Code Length
k = 6   # Code Dimension
r = [4, 3]  # Locality of the code
local_minimum_distance = [5, 5]  # correct one error
sub_group_type = ["add", "mult"]
use_erasure_decoder = False

## Simulation name and directory
sim_name = "GF_" + str(q) + "_" + sub_group_type[0] + "_" + str(r[0]+local_minimum_distance[0]-1) + "_" + sub_group_type[1] + "_" + str(r[1]+local_minimum_distance[1]-1)
if use_erasure_decoder:
    res_dir = './results/' + sim_name + "_erasure"
else:
    res_dir = './results/' + sim_name

res_err_path = res_dir + '/' + sim_name + '_err.csv'

csvfile = open(res_err_path, newline='')

design_distnace = 15
real_distance = 45

succ_count = [0] * real_distance
prob = [0] * real_distance

err_reader = csv.reader(csvfile)

for row in err_reader:
    if row[1] == "remain_err_weight":
        continue
    if (int(row[1]) <= (design_distnace-1 // 2)):
        succ_count[int(row[0])-1] += int(row[2])

print(succ_count)

printed = False
for i in range(0,len(prob)):
    prob[i] = (succ_count[i]/5000)*100
    if ((prob[i] < 99) & (printed == False)):
        print(i-1)
        printed = True

for i in range(0,len(prob)):
    print(i+1, ":", prob[i])