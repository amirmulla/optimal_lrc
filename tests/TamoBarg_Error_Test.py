import sys
sys.path.append('./src')

from sage.all import *
from sage.coding.channel import StaticErrorRateChannel
from TamoBarg import TamoBargCode

# LRC Parameters
q = 16  # Field size
n = 16  # Code dimension
k = 3  # Information dimension
r = 7  # Locality of the code
local_minimum_distance = 2  # Correctable erasures

n_err = 1

print("GF: ", q)
print("Code dim n: ", n)
print("Message dim k: ", k)
print("Locality r: ", r)
print("Local Minimum distance d: ", local_minimum_distance)

# GF
F = GF(q, repr='int')

# Message Space
M = VectorSpace(F, k)
# Choose random message
message = M.random_element()

# Code space
V = VectorSpace(F, n)

C = TamoBargCode(F, n, k, r, local_minimum_distance)

print("Minimum distance d: ", C.minimum_distance())

Dec = C.decoder("ErrorDecoder")
Enc = C.encoder()
c = Enc.encode(message)

Chan = StaticErrorRateChannel(V, n_err)

r = Chan.transmit(c)
e = r - c

print(C)
print("Error              : ", e)
print("Recieved Word      : ", r)
print("Original Codeword  : ", c)

correct_c = Dec.decode_to_code(r)

print("Corrected Codeword : ", correct_c)
print("Correction Successfull: ", correct_c == c)
