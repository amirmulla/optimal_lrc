import sys

sys.path.append('./src')

from sage.coding.channel import ErrorErasureChannel
from sage.all import *
from TamoBerg import *

q = 16  # Field size
n = 16  # Code dimension
k = 8  # Information/message dimension
r = 3  # Locality of the code
local_minimum_distance = 2  # correctable erasures

n_err, n_era = 0, 2

# GF
F = GF(q, repr='int')

# Message Space
M = VectorSpace(F, k)
# Choose random message
message = M.random_element()

# Code space
V = VectorSpace(F, n)

C = TamoBergCode(F, n, k, r, local_minimum_distance)

sub_group, sub_group_type = C.sub_group()
partitions, _ = C.partition()
evalpts = C.evaluation_points()

print("GF: ", q)
print("Code dim n: ", n)
print("Message dim k: ", k)
print("Locality r: ", r)
print("Global Minimum distance d: ", (n - k - ceil(k / r) + 2))
print("Local Minimum distance d: ", local_minimum_distance)
print("Subgroup: ", sub_group)
print("Subgroup Type: ", sub_group_type)
print("Subgroup Size: ", len(sub_group))
print("Partition: ", partitions)
print("Evaluation points: ", evalpts)

Dec = C.decoder()
Enc = C.encoder()
c = Enc.encode(message)

Chan = ErrorErasureChannel(V, n_err, n_era)

r = Chan.transmit(c)

# Perform erasure decoding
try:
    correct_c = Dec.decode_to_code(r)
except:
    correct_c = r[0]

print(C)
print("Erasure Vector: ", r[1])
print("Recieved Word : ", r[0])
print("Codeword      : ", c)
print("Corr Codeword : ", correct_c)
print("Correction Successfull: ", correct_c == c)
