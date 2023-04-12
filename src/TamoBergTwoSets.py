from sage.all import *
from sage.coding.linear_code import AbstractLinearCode
from sage.coding.encoder import Encoder
from sage.coding.decoder import Decoder, DecodingError
from sage.coding.grs_code import GeneralizedReedSolomonCode
import networkx as nx
from networkx.algorithms import *

####################### code ###############################


class TamoBergCodeTwoSets(AbstractLinearCode):
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field, length, dimension, locality, local_minimum_distance=2, sub_group_type="any"):
        super().__init__(base_field, length, "VectorEncoder", "IterativeDecoder")
        if not base_field.is_finite():
            raise ValueError("base_field has to be a finite field")

        if (type(locality) == list):
            self._num_of_sets = len(locality)
        else:
            self._num_of_sets = 1

        self._dimension = dimension
        self._locality = locality
        self._local_minimum_distance = local_minimum_distance
        self._base_mult_group = base_field.list()[1:base_field.cardinality()]
        self._base_mult_group.sort()
        self._mult_sub_groups = self._list_mult_subgroups()
        self._add_sub_groups = self._list_add_subgroups()

        if self._num_of_sets == 1:
            if sub_group_type == "any":
                try:
                    self._mult_sub_groups[locality + local_minimum_distance - 1]
                    self._sub_group = self._mult_sub_groups[locality + local_minimum_distance - 1]
                    self._sub_group_type = "mult"
                except KeyError:
                    try:
                        self._add_sub_groups[locality + local_minimum_distance - 1]
                        self._sub_group = self._add_sub_groups[locality + local_minimum_distance - 1]
                        self._sub_group_type = "add"
                    except KeyError:
                        print("Can't find any subgroup of size: ", locality + local_minimum_distance - 1)

            elif sub_group_type == "mult":
                try:
                    self._mult_sub_groups[locality + local_minimum_distance - 1]
                    self._sub_group = self._mult_sub_groups[locality + local_minimum_distance - 1]
                    self._sub_group_type = "mult"
                except KeyError:
                    print("Can't find mult subgroup of size: ", locality + local_minimum_distance - 1)

            elif sub_group_type == "add":
                try:
                    self._add_sub_groups[locality + local_minimum_distance - 1]
                    self._sub_group = self._add_sub_groups[locality + local_minimum_distance - 1]
                    self._sub_group_type = "add"
                except KeyError:
                    print("Can't find add subgroup of size: ", locality + local_minimum_distance - 1)

            self._parition_size = int(length/(locality + local_minimum_distance - 1))
            self._partition = self._list_sub_group_cosets(sub_group=self._sub_group, sub_group_type=self._sub_group_type)[:self._parition_size]
            self._evaluation_points = flatten(self._partition)
        else:  # Possible combinations of mult and add subgroups
            self._sub_group = []
            self._sub_group_type = sub_group_type
            self._parition_size = []
            self._partition = []
            self._evaluation_points_idx = []
            for i in range(0, self._num_of_sets):
                if (self._sub_group_type[i] == "mult"):
                    try:
                        self._sub_group.append(self._mult_sub_groups[locality[i] + local_minimum_distance[i] - 1])
                    except:
                        print("Can't find mult subgroup of size: ", locality[i] + local_minimum_distance[i] - 1)
                        print("Possible Mult Subgroups: ", self._mult_sub_groups)

                elif (self._sub_group_type[i] == "add"):
                    try:
                        self._sub_group.append(self._add_sub_groups[locality[i] + local_minimum_distance[i] - 1])
                    except:
                        print("Can't find add subgroup of size: ", locality[i] + local_minimum_distance[i] - 1)
                        print("Possible Add Subgroups: ", self._add_sub_groups)

                    # shift group to get the second add subgroup that intersect only at zero
                    if ((i == 1) & (self._sub_group_type[0] == self._sub_group_type[1])):
                        a = base_field.primitive_element()
                        multiplier = (a**log(len(self._sub_group[0]),2))
                        tmp = []
                        for h in self._sub_group[i]:
                            tmp.append(h * multiplier)
                        tmp.sort()
                        self._sub_group[i] = tmp

                self._parition_size.append(int(length /(locality[i] + local_minimum_distance[i] - 1)))
                self._partition.append(self._list_sub_group_cosets(sub_group=self._sub_group[i],sub_group_type=self._sub_group_type[i])[:self._parition_size[i]])

            self._evaluation_points = flatten(self._partition[0])
            self._evaluation_points.sort()

            for i in range(0, self._num_of_sets):
                tmp_a = []
                for coset in self._partition[i]:
                    tmp_b = []
                    for h in coset:
                        tmp_b.append(self._evaluation_points.index(h))
                    tmp_a.append(tmp_b)

                self._evaluation_points_idx.append(tmp_a)

            # Create Graph representation
            self._G, self._G_L, self._G_R = self._create_bipartite_graph()

    def __eq__(self, other):
        return isinstance(other, TamoBergCodeTwoSets) and \
            self.length() == other.length() and \
            self.dimension() == other.dimension() and \
            self.locality() == other.locality()

    def _repr_(self):
        return "[%s, %s, %s] Tamo-Berg Code over GF(%s)" % (self.length(), self.dimension(), self.locality(), self.base_field().cardinality())

    def _list_mult_subgroups(self):
        multiplicative_subgroups = {}
        for j in self._base_mult_group:
            tmp = []
            first_identity_found = 0
            for i in range(0, len(self._base_mult_group)):
                if (j**i == 1 & first_identity_found == 1):
                    break
                tmp.append(j**i)
                first_identity_found = 1

            tmp.sort()
            multiplicative_subgroups[len(tmp)] = tmp

        return multiplicative_subgroups

    def _list_add_subgroups(self):
        additive_subgroups = {}
        F = self.base_field()
        p = F.characteristic()
        for i in range(1, F.degree()):
            tmp = []
            Sub_F = GF(p**i, name="a", repr="int")
            for elm in Sub_F:
                tmp.append(F.from_integer(int(str(elm))))
            tmp.sort()
            additive_subgroups[len(tmp)] = tmp

        return additive_subgroups

    def _list_sub_group_cosets(self, sub_group, sub_group_type):
        F = self.base_field()
        h_cosets = []
        for g in self._base_mult_group:
            tmp = []
            for h in sub_group:
                if (sub_group_type == "mult"):
                    tmp.append(h * g)
                elif (sub_group_type == "add"):
                    tmp.append(h + g)
            tmp.sort()
            h_cosets.append(tmp)

        return list(map(list, set(map(lambda i: tuple(i), h_cosets))))

    def _create_bipartite_graph(self):
        G = nx.Graph()

        # Add graph nodes
        k = 0
        for i in range(0, len(self._partition)):
            for j in range(0, len(self._partition[i])):
                G.add_node(k, partition=i, coset=self._partition[i][j], bipartite=i)
                k += 1

        # Add graph Edges
        for i in G.nodes():
            for e in G.nodes[i]["coset"]:
                for j in G.nodes():
                    if ((e in G.nodes[j]["coset"]) & (i != j)):
                        G.add_edge(i, j, evalpts=e)

        G_L = {n for n, d in G.nodes(data=True) if d["bipartite"] == 0}
        G_R = set(G) - G_L

        return G, G_L, G_R

    def bipartite_graph(self):
        return self._G, self._G_L, self._G_R

    def partition(self):
        return self._partition, self._parition_size

    def sub_group(self):
        return self._sub_group, self._sub_group_type

    def evaluation_points_idx(self):
        return self._evaluation_points_idx

    def evaluation_points(self):
        return self._evaluation_points

    def locality(self):
        return self._locality

    def num_of_sets(self):
        return self._num_of_sets

    def global_minimum_distance(self):
        return (self.length() - self.dimension() - ceil(self.dimension() / self.locality()) + 2)

    def local_minimum_distance(self):
        return self._local_minimum_distance


####################### encoders ###############################


class TamoBergVectorEncoder(Encoder):

    def __init__(self, code):
        self._num_of_sets = code.num_of_sets()
        self._partition, self._partition_size = code.partition()
        self._sub_group, self._sub_group_type = code.sub_group()
        self._base_poly_ring = PolynomialRing(code.base_field(), 'x')
        self._annihilator = self._calc_annihilator(code.evaluation_points(),self._base_poly_ring)
        self._quotient_poly_ring = self._base_poly_ring.quotient(self._annihilator, 'x')
        self._good_poly = []
        self._algebra_basis = []
        self._enc_basis = []
        for i in range(0, self._num_of_sets):
            self._good_poly.append(self._calc_good_poly(self._sub_group[i], self._sub_group_type[i], self._quotient_poly_ring))
            self._algebra_basis.append(self._calc_algebra_basis(self._good_poly[i], self._partition_size[i]))
            self._enc_basis.append(self._calc_enc_basis(self._algebra_basis[i], code.locality()[i], self._quotient_poly_ring))

        self._comb_enc_basis = self._calc_combine_enc_basis(self._enc_basis[0], self._enc_basis[1], self._quotient_poly_ring)

        self._max_dimension = len(self._comb_enc_basis)
        if(code.dimension() > self._max_dimension):
            raise ValueError("code dimension is too big. maximum code dimension: ", self._max_dimension)

        self._comb_enc_basis = self._comb_enc_basis[:code.dimension()]

        super().__init__(code)

    def _repr_(self):
        return "Vector Encoder for %s" % self.code()

    def _calc_annihilator(self, pts, R):
        a = R.gen()
        h = 1
        for p in pts:
            h = h * (a - p)
        return h

    def _calc_good_poly(self, sub_group, sub_group_type, S):
        x = S.gen()
        if (sub_group_type == "mult"):
            good_poly = x**len(sub_group)
        else:
            good_poly = S.one()
            for h in sub_group:
                good_poly = good_poly * (x - h)
        return good_poly

    def _calc_algebra_basis(self, good_poly, partition_size):
        algebra_basis = []
        for i in range(0, partition_size):
            algebra_basis.append(good_poly**i)

        return algebra_basis

    def _calc_enc_basis(self, algebra_basis, r, S):
        x = S.gen()
        enc_basis = []
        for i in range(0, r):
            for p in algebra_basis:
                enc_basis.append(p * (x**i))

        enc_basis.sort()

        return enc_basis

    def _calc_combine_enc_basis(self, a_enc_basis, b_enc_basis, S):
        F = S.base_ring()

        U1 = matrix(F, [s.list() for s in a_enc_basis]).transpose()
        U2 = matrix(F, [s.list() for s in b_enc_basis]).transpose()
        U2 = -U2

        A = copy(U1)
        A = A.augment(U2)
        M = A.right_kernel_matrix()
        E = M.matrix_from_columns(list(range(0, U1.ncols())))
        T = (E * U1.transpose())

        x = S.gen()
        comb_enc_basis = []
        for row in T:
            tmp = S.zero()
            for i in range(0,len(row)):
                tmp += row[i] * (x**i)
            if (tmp != S.zero()):
                comb_enc_basis.append(tmp)

        comb_enc_basis.sort()
        return comb_enc_basis

    def enc_basis_ext(self):
        return self._enc_basis

    def enc_algebra_basis(self):
        return self._algebra_basis

    def enc_good_poly(self):
        return self._good_poly

    def enc_basis(self):
        return self._comb_enc_basis

    def enc_poly(self):
        return self._enc_poly

    @cached_method
    def generator_matrix(self):
        C = self.code()
        alphas = C.evaluation_points()
        g = matrix(C.base_field(), C.dimension(), C.length(),
                   lambda i, j: self._comb_enc_basis[i].lift()(alphas[j]))
        g.set_immutable()
        return g


######################### decoders #################################

###### TamoBergIterariveDecoder #######
class TamoBergIterariveDecoder(Decoder):

    def __init__(self, code, max_num_of_itr=2):
        input_space = cartesian_product([code.ambient_space(), VectorSpace(GF(2), code.ambient_space().dimension())])
        super().__init__(code, input_space, "VectorEncoder")
        self._max_num_of_itr = max_num_of_itr

    def _repr_(self):
        return "Iterative Decoder for %s" % self.code()

    def decode_to_code(self, r):
        C = self.code()
        F = C.base_field()
        partitions, _ = C.partition()
        evalpts_idx = C.evaluation_points_idx()

        i = 0
        partitions_local_decoders = []
        partitions_local_codes = []
        for partition in partitions:
            code_tmp = []
            decoder_tmp = []
            for coset in partition:
                local_grs = GeneralizedReedSolomonCode(coset, C.locality()[i])
                code_tmp.append(local_grs)
                decoder_tmp.append(local_grs.decoder("Gao"))

            partitions_local_decoders.append(decoder_tmp)
            partitions_local_codes.append(code_tmp)
            i += 1

        uncorr_err = True
        num_itr = 0
        while ((uncorr_err is True) & (num_itr < self._max_num_of_itr)):
            uncorrectable_err = [False, False]
            for i in range(0,len(partitions)):
                partition = partitions[i]
                for j in range(0,len(partition)):
                    coset = partition[j]
                    r_list = []
                    for k in range(0,len(coset)):
                        r_list.append(r[evalpts_idx[i][j][k]])

                    try:
                        corr_c = partitions_local_decoders[i][j].decode_to_code(vector(F, r_list))
                    except: # Decoding fails
                        uncorrectable_err[i] = True
                        corr_c = vector(F, r_list)

                    r_list = list(corr_c)

                    # Fix error in r
                    for l in range(0,len(coset)):
                        r[evalpts_idx[i][j][l]] = r_list[l]

            if (uncorrectable_err[0] == False & uncorrectable_err[1] == False):
                uncorr_err = False
            num_itr += 1

        return r, num_itr

###### Error-and-Erasure Iterarive Decoder #######
class TamoBergIterariveEEDecoder(Decoder):

    def __init__(self, code, max_num_of_itr=2):
        input_space = cartesian_product([code.ambient_space(), VectorSpace(GF(2), code.ambient_space().dimension())])
        super().__init__(code, input_space, "VectorEncoder")
        self._max_num_of_itr = max_num_of_itr

    def _repr_(self):
        return "Iterative Decoder Error-and-Erasure for %s" % self.code()

    def decode_to_code(self, r, orig_c):
        print_log = False
        C = self.code()
        F = C.base_field()
        partitions, _ = C.partition()
        evalpts_idx = C.evaluation_points_idx()

        i = 0
        partitions_local_decoders = []
        partitions_local_codes = []
        for partition in partitions:
            code_tmp = []
            decoder_tmp = []
            for coset in partition:
                local_grs = GeneralizedReedSolomonCode(coset, C.locality()[i])
                code_tmp.append(local_grs)
                decoder_tmp.append(local_grs.decoder("ErrorErasure"))

            partitions_local_decoders.append(decoder_tmp)
            partitions_local_codes.append(code_tmp)
            i += 1


        erasure_start_itr = 5

        erasure_vector = zero_vector(GF(2), len(r))
        uncorr_err = True
        num_itr = 0
        while ((uncorr_err is True) & (num_itr < self._max_num_of_itr)):
            uncorrectable_err = [False, False]
            due_status = []
            for i in range(0,len(partitions)):
                partition = partitions[i]
                due_tmp = []
                for j in range(0,len(partition)):
                    coset = partition[j]
                    r_list = []
                    orig_c_list = []
                    erasure_vector_list = []
                    for k in range(0,len(coset)):
                        r_list.append(r[evalpts_idx[i][j][k]])
                        orig_c_list.append(orig_c[evalpts_idx[i][j][k]])
                        erasure_vector_list.append(erasure_vector[evalpts_idx[i][j][k]])


                    word_and_erasure_vector = vector(F, r_list), vector(GF(2), erasure_vector_list)

                    erasure_vec = vector(GF(2), erasure_vector_list)

                    try:
                        corr_c = partitions_local_decoders[i][j].decode_to_code(word_and_erasure_vector)
                        due_tmp.append(False)
                    except: # Decoding fails
                        due_tmp.append(True)
                        uncorrectable_err[i] = True
                        corr_c = vector(F, r_list)

                    if ((num_itr >= erasure_start_itr) & (print_log == True)):
                        print("#######################################################")
                        print(j,":",word_and_erasure_vector)
                        print("Orig_c:" ,orig_c_list)
                        print("Error:", vector(F, r_list) - vector(F, orig_c_list))
                        print("corr_c:", corr_c)
                        print("DUE Status:", due_tmp[j])
                        print(partitions_local_decoders[i][j])
                        print("#######################################################")

                    if erasure_vec.hamming_weight() == partitions_local_codes[i][j].minimum_distance()-1:
                        pts = []
                        err = []
                        for k in range(0, len(coset)):
                            if (erasure_vec[k] == 0):
                                pts.append((coset[k], r_list[k]))
                            else:
                                err.append(k)

                        R = PolynomialRing(C.base_field(), 'x')
                        d = R.lagrange_polynomial(pts)

                        corr_c = vector(F, r_list)
                        for e in err:
                            corr_c[e] = d(coset[e])

                    r_list = list(corr_c)

                    # Fix error in r
                    for l in range(0,len(coset)):
                        r[evalpts_idx[i][j][l]] = r_list[l]

                due_status.append(due_tmp)

            if (uncorrectable_err[0] == False & uncorrectable_err[1] == False):
                uncorr_err = False

            if print_log:
                print(num_itr, ":", due_status)

            num_itr += 1

            if (num_itr >= erasure_start_itr):
                for s in range(0,len(due_status)):
                    p = partitions[s]
                    for t in range(0,len(p)):
                        for m in range(0, len(p[t])):
                            if((due_status[s][t] == True) & (s == 0)):
                                erasure_vector[evalpts_idx[s][t][m]] = 1
                            elif(due_status[s][t] == False):
                                erasure_vector[evalpts_idx[s][t][m]] = 0

                if print_log:
                    print("Erasure Vector:", erasure_vector)
                    print("Error         :", r-orig_c)

        return r, num_itr




####################### registration ###############################
TamoBergCodeTwoSets._registered_encoders["VectorEncoder"] = TamoBergVectorEncoder
TamoBergCodeTwoSets._registered_decoders["IterativeDecoder"] = TamoBergIterariveDecoder
TamoBergCodeTwoSets._registered_decoders["IterativeErasureErrorDecoder"] = TamoBergIterariveEEDecoder
