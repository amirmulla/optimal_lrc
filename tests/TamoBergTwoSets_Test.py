import sys

sys.path.append('./src')

from sage.all import *
from TamoBergTwoSets import *
from sage.coding.channel import StaticErrorRateChannel

q = 64  # Field size
n = 63  # Code dimension
k = 35 # Information/message dimension
r = [5, 7]  # Locality of the code
local_minimum_distance = [3, 3]  # correct one error
sub_group_type = ["mult", "mult"]
n_err = 5
max_num_of_itr = 10

print("GF: ", q)
print("Code dim n: ", n)
print("Message dim k: ", k)
print("Locality r1: ", r[0])
print("Locality r2: ", r[1])
print("Local Minimum distance d1: ", local_minimum_distance[0])
print("Local Minimum distance d2: ", local_minimum_distance[1])

# GF
F = GF(q, repr='int')

# Code space
V = VectorSpace(F, n)

# Message Space
M = VectorSpace(F, k)
# Choose random message
message = M.random_element()

C = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type)
print(C)
Dec = C.decoder("IterativeDecoder", max_num_of_itr)
Enc = C.encoder()
c = Enc.encode(message)

partitions, _ = C.partition()
evalpts = C.evaluation_points()

Chan = StaticErrorRateChannel(V, n_err)

r = Chan.transmit(c)
e = r - c

print("Error              : ", e)
print("Recieved Word      : ", r)
print("Original Codeword  : ", c)

correct_c = Dec.decode_to_code(r)

print("Corrected Codeword : ", correct_c)
print("Correction Successfull: ", correct_c[0] == c)
