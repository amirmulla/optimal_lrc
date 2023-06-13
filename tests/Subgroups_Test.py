from sage.modules.free_module import VectorSpace
from sage.all import *
from itertools import combinations
import sys

sys.path.append('./src')


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
            subgroup.append(F.from_integer(int(binary_string, F.characteristic())))

        subgroup.sort()
        subgroups.append(subgroup)

    return subgroups


q = 81  # Field size
print("GF: ", q)
# GF
F = GF(q, name='a', repr='int')

additive_subgroups = find_additive_subgroups(F, 9)

base_mult_group = F.list()[1:]
multiplicative_subgroups = {}
for j in base_mult_group:
    tmp = []
    first_identity_found = 0
    for i in range(0, len(base_mult_group)):
        if (j**i == 1 & first_identity_found == 1):
            break
        tmp.append(j**i)
        first_identity_found = 1

    tmp.sort()
    multiplicative_subgroups[len(tmp)] = tmp


def list_sub_group_cosets(sub_group, sub_group_type):
    h_cosets = []
    for g in base_mult_group:
        tmp = []
        for h in sub_group:
            if (sub_group_type == "mult"):
                tmp.append(h * g)
            elif (sub_group_type == "add"):
                tmp.append(h + g)
        tmp.sort()
        h_cosets.append(tmp)

    return list(map(list, set(map(lambda i: tuple(i), h_cosets))))


print(multiplicative_subgroups)
for additive_subgroup in additive_subgroups:
    sub_group_0 = multiplicative_subgroups[4]
    sub_group_1 = additive_subgroup  # additive_subgroups[6]

    #print(sub_group_0)
    #print(sub_group_1)

    #print("###########################")

    m_costs = list_sub_group_cosets(sub_group_0, "mult")
    a_costs = list_sub_group_cosets(sub_group_1, "add")

    #for x in m_costs:
    #    print(x)

    #print("###########################")

    #for y in a_costs:
    #    print(y)


    def intersection(lst1, lst2):
        lst3 = [value for value in lst1 if value in lst2]
        return lst3

    chk = 0
    for x in m_costs:
        for y in a_costs:
            tmp = intersection(x, y)
            if len(tmp) > 1:
                chk = 1

    if chk == 0:
        print(sub_group_0)
        print(sub_group_1)
        print("###########################")
