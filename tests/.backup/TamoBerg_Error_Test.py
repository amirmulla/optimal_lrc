import sys
sys.path.append('./src')

from sage.all import *
from sage.coding.channel import StaticErrorRateChannel
from TamoBerg import TamoBergCode

# LRC Parameters
q = 243  # Field size
n = 16  # Code dimension
k = 8  # Information dimension
r = 2  # Locality of the code
local_minimum_distance = 3  # Correctable erasures

n_err = 1


print(factor(q))
print(factor(q-1))
exit()
print("GF: ", q)
print("Code dim n: ", n)
print("Message dim k: ", k)
print("Locality r: ", r)
print("Local Minimum distance d: ", local_minimum_distance)

# GF
F = GF(q, name='alpha', repr='log')

# Message Space
M = VectorSpace(F, k)
# Choose random message
message = M.random_element()

# Code space
V = VectorSpace(F, n)

C = TamoBergCode(F, n, k, r, local_minimum_distance, sub_group_type="add")

print("Minimum distance d: ", C.minimum_distance())

Dec = C.decoder("ErrorDecoder")
Enc = C.encoder()
c = Enc.encode(message)

Chan = StaticErrorRateChannel(V, n_err)

r = Chan.transmit(c)
e = r - c

print("Partition          : ", C.partition())
print("Subgroup           : ", C.sub_group())
print("Good Poly          : ", Enc.enc_good_poly())
print("Encode Algebra     : ", Enc.enc_algebra_basis())
print("Encode Basis       : ", Enc.enc_basis())
print("Message            : ", message)
print("Evaluation points  : ", C.evaluation_points())
print("Error              : ", e)
print("Recieved Word      : ", r)
print("Original Codeword  : ", c)


print("Partition          : ", latex(C.partition()))
print("Subgroup           : ", latex(C.sub_group()))
print("Good Poly          : ", latex(Enc.enc_good_poly()))
print("Encode Algebra     : ", latex(Enc.enc_algebra_basis()))
print("Encode Basis       : ", latex(Enc.enc_basis()))
print("Message            : ", latex(message))
print("Evaluation points  : ", latex(C.evaluation_points()))
print("Error              : ", latex(e))
print("Recieved Word      : ", latex(r))
print("Original Codeword  : ", latex(c))


correct_c = Dec.decode_to_code(r)

print("Corrected Codeword : ", correct_c)
print("Correction Successfull: ", correct_c == c)
