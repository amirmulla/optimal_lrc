import sys
import time
sys.path.append('./src')

from TamoBergTwoSets import TamoBergCodeTwoSets
from sage.coding.channel import StaticErrorRateChannel
from sage.all import *

#############


def is_closed_under_addition(subgroup):
    """
    Check if a subgroup is closed under addition.

    Args:
        subgroup: A set representing the subgroup.

    Returns:
        True if the subgroup is closed under addition, False otherwise.
    """
    for a in subgroup:
        for b in subgroup:
            if (a + b) not in subgroup:
                return False

    return True


def find_additive_subgroups_expo(F, n, num_of_subgroups=1000):
    """
    Find all subgroups of size n in a finite field.

    Args:
        F: A SageMath finite field object.
        n: The size of the subgroups.

    Returns:
        A list of sets representing the subgroups.
    """
    subgroups = []

    # Generate all possible elements in the finite field
    elements = F.list()

    # Generate all possible combinations of elements of size n
    combinations_list = list(combinations(elements, n))

    # Check if each combination forms a subgroup
    for combo in combinations_list:
        subgroup = list(combo)

        # Check subgroup closure under addition
        if is_closed_under_addition(subgroup):
            subgroup.sort()
            subgroups.append(subgroup)

        if len(subgroups) > num_of_subgroups:
            break

    return subgroups


def find_subspaces(n, m, F):
    """
    Find all subspaces of size n in a finite field of order 2^m.

    Parameters:
        n (int): The size of the subspaces.
        m (int): The exponent of the finite field.
        F: A SageMath finite field object.

    Returns:
        list: A list of subspaces, where each subspace is represented by its basis.
    """
    # Create the finite field of order 2^m
    k = F.characteristic()
    Sub_F = GF(k, 'a', repr='int')

    # Create the vector space over the finite field
    V = VectorSpace(Sub_F, m)

    # Find all subspaces of size n
    subspaces = list(V.subspaces(n))

    return subspaces


def find_additive_subgroups(F, n):
    """
    Find all subgroups of size n in a finite field.

    Args:
        F: A SageMath finite field object.
        n: The size of the subgroups.

    Returns:
        A list of sets representing the subgroups.
    """
    subgroups = []

    subspaces = find_subspaces(log(n, F.characteristic()), F.degree(), F)

    for subspace in subspaces:
        subgroup = []
        for vec in subspace.list():
            binary_string = ''.join(str(i) for i in vec)
            subgroup.append(F.from_integer(
                int(binary_string, F.characteristic())))

        subgroup.sort()
        subgroups.append(subgroup)

    return subgroups

##############

# set_random_seed(100)

q = 64  # Field size
n = 64  # Code dimension
k = 6  # Information/message dimension
r = [4, 3]  # Locality of the code
local_minimum_distance = [5, 6]  # correct one error
sub_group_type = ["add", "add"]
n_err = 1
max_num_of_itr = 10

# GF
F = GF(q, name="a", repr='int')

additive_subgroups = find_additive_subgroups(F, 8)
add_subgroup = additive_subgroups[5]

subgroup = [additive_subgroups[3], additive_subgroups[100]]

print(subgroup)
# Code space
V = VectorSpace(F, n)

# Message Space
M = VectorSpace(F, k)
# Choose random message
message = M.random_element()
# message = zero_vector(F,k)

C = TamoBergCodeTwoSets(F, n, k, r, local_minimum_distance, sub_group_type, subgroup=subgroup)

sub_group, sub_group_type = C.sub_group()
partitions, _ = C.partition()
evalpts = C.evaluation_points()

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
