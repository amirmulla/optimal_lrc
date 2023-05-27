import sys

sys.path.append('./src')

from sage.all import *

def unique(list1):
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return unique_list

# LRC Parameters
q = 16  # Field size

print("GF: ", q)

# GF
F = GF(q, name='a', repr='int')

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

print(multiplicative_subgroups)

additive_subgroups = {}
p = F.characteristic()
for i in range(1, F.degree()):
    tmp = []
    Sub_F = GF(p**i, name="a", repr="int")
    for elm in Sub_F:
        tmp.append(F.from_integer(int(str(elm))))

    tmp.sort()
    additive_subgroups[len(tmp)] = tmp

print(additive_subgroups)


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

sub_group_0 = multiplicative_subgroups[3]
sub_group_1 = additive_subgroups[4]

a = F.primitive_element()
p = F.characteristic()
multiplier = (a**log(len(sub_group_1),p))
tmp = []
for h in sub_group_1:
    tmp.append(h * multiplier)

tmp.sort()
sub_group_1 = tmp

print(sub_group_0)
print(sub_group_1)

print("###########################")

m_costs = list_sub_group_cosets(sub_group_0, "mult")
a_costs = list_sub_group_cosets(sub_group_1, "add")

for x in m_costs:
    print(x)

print("###########################")

for y in a_costs:
    print(y)


