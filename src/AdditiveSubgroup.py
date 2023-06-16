from itertools import combinations
from sage.all import *

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
