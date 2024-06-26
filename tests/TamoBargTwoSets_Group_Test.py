from sage.all import *
from itertools import combinations
import sys
import time

sys.path.append('./src')
import networkx as nx
import matplotlib.pyplot as plt
from networkx.algorithms import *
from TamoBargTwoSets import TamoBargCodeTwoSets
from sage.coding.channel import StaticErrorRateChannel
from AdditiveSubgroup import *

# set_random_seed(100)

q = 64  # Field size
n = 64  # Code dimension
k = 12  # Information/message dimension
r = [2, 5]  # Locality of the code
local_minimum_distance = [3, 3]  # correct two errors
sub_group_type = ["add", "mult"]
n_err = 1
max_num_of_itr = 10

# GF
F = GF(q, name="a", repr='int')

# Specify Additive subgroup
if sub_group_type[0] == "add":
    sub_group_size = r[0] + local_minimum_distance[0] -1
    additive_subgroups = find_additive_subgroups(F, sub_group_size)

#add_subgroup = [None, None]
#C = TamoBargCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type, shift_add=False, subgroup=add_subgroup)

min_dd = n
max_dd = 0

min_k = n
max_k = 0

i = 0
for g in additive_subgroups:
    add_subgroup = [g, None]
    C = TamoBargCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type, shift_add=False, subgroup=add_subgroup)

    if C.max_dimension() > max_k:
        max_k = C.max_dimension()
        print("max_k", max_k)
        print("max_k idx", i)

    if C.max_dimension() < min_k:
        min_k = C.max_dimension()
        print("min_k", min_k)
        print("min_k idx", i)

    if C.design_distance() > max_dd:
        max_dd = C.design_distance()
        print("max_dd", max_dd)
        print("max_dd idx", i)

    if C.design_distance() < min_dd:
        min_dd = C.design_distance()
        print("min_dd", min_dd)
        print("min_dd idx", i)

    i += 1

# Code space
V = VectorSpace(F, n)

# Message Space
M = VectorSpace(F, k)
# Choose random message
message = M.random_element()
# message = zero_vector(F,k)

sub_group, sub_group_type = C.sub_group()
partitions, _ = C.partition()
evalpts = C.evaluation_points()

for p in partitions:
    i = 1
    for coset in p:
        print(i, ": ",coset)
        i += 1
    print("########")

G, G_L, G_R = C.bipartite_graph()
m, n = len(G_L), len(G_R)
pos = dict()
pos.update((i, (i - m / 2, 1)) for i in range(m))
pos.update((i, (i - m - n / 2, 0)) for i in range(m, m + n))
nx.draw(G, with_labels=True, pos=pos, node_size=300, width=0.4)
#plt.savefig('graph.png')
plt.show()

#print(G.nodes[7])
#print(G.nodes[11])

print("GF: ", q)
print("Code dim n: ", n)
print("Message dim k: ", k)
print("max k: ", C.max_dimension())
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

# åDec = C.decoder("IterativeErasureErrorDecoder", max_num_of_itr)
Dec = C.decoder("TwoStepsErasureErrorDecoder", max_num_of_itr)
Enc = C.encoder()
c = Enc.encode(message)

#print("######## Encoder Debug ##########")
#enc_basis = Enc.enc_basis_ext()
#enc_algebra_basis = Enc.enc_algebra_basis()
#enc_good_poly = Enc.enc_good_poly()
#enc_comb_enc_basis = Enc.enc_basis()
#print("enc_good_poly: ", enc_good_poly[0])
#print("enc_algebra_basis: ", enc_algebra_basis[0])
#print("enc_basis: ", enc_basis[0])
#print("---------------------------------------")
#print("enc_good_poly: ", enc_good_poly[1])
#print("enc_algebra_basis: ", enc_algebra_basis[1])
#print("enc_basis: ", enc_basis[1])
#print("---------------------------------------")
#print("enc_comb_enc_basis: ", enc_comb_enc_basis)
Chan = StaticErrorRateChannel(V, n_err)
num_of_mrs = 1

print_log = True

st = time.time()
for i in range(0, num_of_mrs):
    r = Chan.transmit(c)
    e = r - c
    correct_c, num_of_itr, num_of_rem_era = Dec.decode_to_code(r)
    if print_log:
        print("Error              : ", e)
        print("Recieved Word      : ", r)
        print("Corrected Codeword : ", correct_c)
        print("Correction Successfull: ", correct_c == c)
        print("Num of Iterations: ", num_of_itr)
        print("Num of Remaining Erasures: ", num_of_rem_era)

et = time.time()
elapsed_time = (et - st) * (10 ** 3)

print('Average Decoding Function Execution time:',
      (elapsed_time/num_of_mrs), 'msec')
