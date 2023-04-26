import sys

sys.path.append('./src')

from sage.all import *
from TamoBergTwoSets import *
from sage.coding.channel import StaticErrorRateChannel

#set_random_seed(100)

q = 31  # Field size
n = 30  # Code dimension
k = 12  # Information/message dimension
r = [3, 4]  # Locality of the code
local_minimum_distance = [3, 3]  # correct one error
sub_group_type = ["mult", "mult"]
n_err = 8
max_num_of_itr = 10

# GF
F = GF(q, name="a", repr='int')

# Code space
V = VectorSpace(F, n)

# Message Space
M = VectorSpace(F, k)
# Choose random message
message = M.random_element()
#message = zero_vector(F,k)

C = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type)
print(C)
sub_group, sub_group_type = C.sub_group()
partitions, _ = C.partition()
evalpts = C.evaluation_points()

print("GF: ", q)
print("Code dim n: ", n)
print("Message dim k: ", k)
print("Locality r1: ", r[0])
print("Locality r2: ", r[1])
print("Local Minimum distance d1: ", local_minimum_distance[0])
print("Local Minimum distance d2: ", local_minimum_distance[1])
print("First Subgroup: ", sub_group[0])
print("First Subgroup Type: ", sub_group_type[0])
print("First Subgroup Size: ", len(sub_group[0]))
print("Second Subgroup: ", sub_group[1])
print("Second Subgroup Type: ", sub_group_type[1])
print("Second Subgroup Size: ", len(sub_group[1]))

Dec = C.decoder("IterativeErasureErrorDecoder", max_num_of_itr)
Enc = C.encoder()
c = Enc.encode(message)

print("######## Encoder Debug ##########")

enc_basis = Enc.enc_basis_ext()
enc_algebra_basis = Enc.enc_algebra_basis()
enc_good_poly = Enc.enc_good_poly()
enc_comb_enc_basis = Enc.enc_basis()

print("enc_good_poly: ", enc_good_poly[0])
print("enc_algebra_basis: ", enc_algebra_basis[0])
print("enc_basis: ", enc_basis[0])
print("---------------------------------------")
print("enc_good_poly: ", enc_good_poly[1])
print("enc_algebra_basis: ", enc_algebra_basis[1])
print("enc_basis: ", enc_basis[1])
print("---------------------------------------")
print("enc_comb_enc_basis: ", enc_comb_enc_basis)

Chan = StaticErrorRateChannel(V, n_err)

r = Chan.transmit(c)
e = r - c

print("Error              : ", e)
print("Recieved Word      : ", r)
print("Original Codeword  : ", c)

correct_c = Dec.decode_to_code(r,c)

print("Corrected Codeword : ", correct_c)
print("Correction Successfull: ", correct_c[0] == c)
print("evalpts: ", evalpts)
print("Error: ", correct_c[0] - c)