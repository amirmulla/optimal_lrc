import sys
import time
sys.path.append('./src')

from sage.all import *
from sage.coding.channel import StaticErrorRateChannel
from TamoBargTwoSets import TamoBargCodeTwoSets
from AdditiveSubgroup import *

q = 81  # Field size
n = 81  # Code dimension
k = 9 # Information/message dimension
r = [3, 3]  # Locality of the code
local_minimum_distance = [7, 7]  # correct one error
sub_group_type = ["add", "add"]
n_err = 1
max_num_of_itr = 10

# GF
F = GF(q, name="alpha", repr='int')

# Code space
V = VectorSpace(F, n)

# Message Space
M = VectorSpace(F, k)
# Choose random message
message = M.random_element()

#C = TamoBargCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type)

# Specify Additive subgroup
#if sub_group_type[0] == "add":
#    sub_group_size = r[0] + local_minimum_distance[0] -1
#    additive_subgroups = find_additive_subgroups(F, sub_group_size)

#add_subgroup = [additive_subgroups[0], None]
C = TamoBargCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type, shift_add=False)#, subgroup=add_subgroup)

sub_group, sub_group_type = C.sub_group()
#partitions, _ = C.partition()

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
exit()

# Dec = C.decoder("IterativeDecoder", max_num_of_itr)
Dec = C.decoder("TwoStepsDecoder", max_num_of_itr)
#Dec = C.decoder("GlobalDecoder")
Enc = C.encoder()
c = Enc.encode(message)

Chan = StaticErrorRateChannel(V, n_err)

r = Chan.transmit(c)
e = r - c

print("Partition          : ", C.partition())
print("Subgroup           : ", C.sub_group())
print("Good Poly          : ", Enc.enc_good_poly())
print("Encode Algebra     : ", Enc.enc_algebra_basis())
print("Encode Basis       : ", Enc.enc_basis_ext())
print("Combined Encode B  : ", Enc.enc_basis())
print("Message            : ", message)
print("Evaluation points  : ", C.evaluation_points())
print("Error              : ", e)
print("Recieved Word      : ", r)
print("Original Codeword  : ", c)


print("Partition          : ", latex(C.partition()))
print("Subgroup           : ", latex(C.sub_group()))
print("Good Poly          : ", latex(Enc.enc_good_poly()))
print("Encode Algebra     : ", latex(Enc.enc_algebra_basis()))
print("Encode Basis       : ", latex(Enc.enc_basis_ext()))
print("Comb Encode Basis  : ", latex(Enc.enc_basis()))
print("Message            : ", latex(message))
print("Evaluation points  : ", latex(C.evaluation_points()))
print("Error              : ", latex(e))
print("Recieved Word      : ", latex(r))
print("Original Codeword  : ", latex(c))


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


partitions, _ = C.partition()
evalpts = C.evaluation_points()

Chan = StaticErrorRateChannel(V, n_err)

avg = 0
num_of_mrs = 10
print_log = True

st = time.time()
for i in range(0, num_of_mrs):
    r = Chan.transmit(c)
    e = r - c
    correct_c, _ = Dec.decode_to_code(r)

    if print_log:
        print("Error              : ", e)
        print("Recieved Word      : ", r)
        print("Corrected Codeword : ", correct_c)
        print("Correction Successfull: ", correct_c == c)

et = time.time()
elapsed_time = (et - st) * (10 ** 3)
print('Average Decoding Function Execution time:', (elapsed_time/num_of_mrs), 'msec')
