from sage.all import *
from sage.coding.linear_code import AbstractLinearCode
from sage.coding.encoder import Encoder
from sage.coding.decoder import Decoder, DecodingError
from sage.coding.grs_code import GeneralizedReedSolomonCode
import networkx as nx
from networkx.algorithms import *

####################### code ###############################


class TamoBargCodeTwoSets(AbstractLinearCode):
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field, length, dimension, locality, local_minimum_distance=2, sub_group_type="any", shift_add=False, subgroup=[None, None]):
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
                    self._mult_sub_groups[locality +
                                          local_minimum_distance - 1]
                    self._sub_group = self._mult_sub_groups[locality +
                                                            local_minimum_distance - 1]
                    self._sub_group_type = "mult"
                except KeyError:
                    try:
                        self._add_sub_groups[locality +
                                             local_minimum_distance - 1]
                        self._sub_group = self._add_sub_groups[locality +
                                                               local_minimum_distance - 1]
                        self._sub_group_type = "add"
                    except KeyError:
                        print("Can't find any subgroup of size: ",
                              locality + local_minimum_distance - 1)

            elif sub_group_type == "mult":
                try:
                    self._mult_sub_groups[locality +
                                          local_minimum_distance - 1]
                    self._sub_group = self._mult_sub_groups[locality +
                                                            local_minimum_distance - 1]
                    self._sub_group_type = "mult"
                except KeyError:
                    print("Can't find mult subgroup of size: ",
                          locality + local_minimum_distance - 1)

            elif sub_group_type == "add":
                try:
                    self._add_sub_groups[locality + local_minimum_distance - 1]
                    self._sub_group = self._add_sub_groups[locality +
                                                           local_minimum_distance - 1]
                    self._sub_group_type = "add"
                except KeyError:
                    print("Can't find add subgroup of size: ",
                          locality + local_minimum_distance - 1)

            self._parition_size = int(
                length/(locality + local_minimum_distance - 1))
            self._partition = self._list_sub_group_cosets(
                sub_group=self._sub_group, sub_group_type=self._sub_group_type)[:self._parition_size]
            self._evaluation_points = flatten(self._partition)
        else:  # Possible combinations of mult and add subgroups
            self._sub_group = []
            self._sub_group_type = sub_group_type
            self._parition_size = []
            self._partition = []
            self._evaluation_points_idx = []
            for i in range(0, self._num_of_sets):
                if(subgroup[i] == None):
                    if (self._sub_group_type[i] == "mult"):
                        try:
                            self._sub_group.append(self._mult_sub_groups[locality[i] + local_minimum_distance[i] - 1])
                        except:
                            print("Can't find mult subgroup of size: ",locality[i] + local_minimum_distance[i] - 1)
                            print("Possible Mult Subgroups: ",self._mult_sub_groups)

                    elif (self._sub_group_type[i] == "add"):
                        try:
                            self._sub_group.append(self._add_sub_groups[locality[i] + local_minimum_distance[i] - 1])
                        except:
                            print("Can't find add subgroup of size: ",locality[i] + local_minimum_distance[i] - 1)
                            print("Possible Add Subgroups: ", self._add_sub_groups)

                        # shift group to get the second add subgroup that intersect only at zero
                        if ((i == 1) & (self._sub_group_type[0] == self._sub_group_type[1])) or shift_add:
                            a = base_field.primitive_element()
                            p = base_field.characteristic()
                            multiplier = (a**log(len(self._sub_group[0]), p))
                            tmp = []
                            for h in self._sub_group[i]:
                                tmp.append(h * multiplier)
                            tmp.sort()
                            self._sub_group[i] = tmp
                else:
                    self._sub_group.append(subgroup[i])

                self._parition_size.append(int(length / (locality[i] + local_minimum_distance[i] - 1)))
                self._partition.append(self._list_sub_group_cosets(sub_group=self._sub_group[i], sub_group_type=self._sub_group_type[i])[:self._parition_size[i]])

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

            if (self.dimension() > self.max_dimension()):
                raise ValueError(
                    "code dimension is too big. maximum code dimension: ", self.max_dimension())

    def __eq__(self, other):
        return isinstance(other, TamoBargCodeTwoSets) and \
            self.length() == other.length() and \
            self.dimension() == other.dimension() and \
            self.locality() == other.locality()

    def _repr_(self):
        return "[%s, %s, %s] Tamo-Barg Code over GF(%s)" % (self.length(), self.dimension(), self.locality(), self.base_field().cardinality())

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
        G = nx.MultiGraph()

        # Add graph nodes
        k = 0
        for i in range(0, len(self._partition)):
            for j in range(0, len(self._partition[i])):
                G.add_node(k, partition=i, coset=self._partition[i][j], bipartite=i)
                k += 1

        # Add graph Edges
        evalpts_done = []
        for i in G.nodes():
            for e in G.nodes[i]["coset"]:
                if e not in evalpts_done:
                    for j in G.nodes():
                        if ((e in G.nodes[j]["coset"]) & (i != j)):
                            G.add_edge(i, j, evalpts=e)
                            evalpts_done.append(e)

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

    def design_distance(self):
        Enc = self.encoder()
        max_degree = 0
        for poly in Enc.enc_basis():
            if (poly.lift().degree() > max_degree):
                max_degree = poly.lift().degree()
        return (self.length() - max_degree)

    def max_dimension(self):
        Enc = self.encoder()
        return Enc.max_dimension()

    def max_rate(self):
        return self.max_dimension() / self.length()

    def local_minimum_distance(self):
        return self._local_minimum_distance

    def minimum_distance_py(self):
        G = self.generator_matrix()
        M = VectorSpace(self.base_field(), self.dimension())
        zero_vector = vector(self.base_field(), [
                             self.base_field().zero()] * self.length())
        d = self.length()
        for m in M:
            c = m*G
            if (c != zero_vector):
                if (c.hamming_weight() < d):
                    d = c.hamming_weight()

        return d


####################### encoders ###############################


class TamoBargVectorEncoder(Encoder):

    def __init__(self, code):
        self._num_of_sets = code.num_of_sets()
        self._partition, self._partition_size = code.partition()
        self._sub_group, self._sub_group_type = code.sub_group()
        self._base_poly_ring = PolynomialRing(code.base_field(), 'x')
        self._annihilator = self._calc_annihilator(
            code.evaluation_points(), self._base_poly_ring)
        self._quotient_poly_ring = self._base_poly_ring.quotient(
            self._annihilator, 'x')
        self._good_poly = []
        self._algebra_basis = []
        self._enc_basis = []
        for i in range(0, self._num_of_sets):
            self._good_poly.append(self._calc_good_poly(
                self._sub_group[i], self._sub_group_type[i], self._quotient_poly_ring))
            self._algebra_basis.append(self._calc_algebra_basis(
                self._good_poly[i], self._partition_size[i]))
            self._enc_basis.append(self._calc_enc_basis(
                self._algebra_basis[i], code.locality()[i], self._quotient_poly_ring))

        self._comb_enc_basis = self._calc_combine_enc_basis(
            self._enc_basis[0], self._enc_basis[1], self._quotient_poly_ring)

        self._max_dimension = len(self._comb_enc_basis)
        if (code.dimension() > self._max_dimension):
            raise ValueError(
                "code dimension is too big. maximum code dimension: ", self._max_dimension)

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
            for i in range(0, len(row)):
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

    def max_dimension(self):
        return self._max_dimension

    @cached_method
    def generator_matrix(self):
        C = self.code()
        alphas = C.evaluation_points()
        g = matrix(C.base_field(), C.dimension(), C.length(),
                   lambda i, j: self._comb_enc_basis[i].lift()(alphas[j]))
        g.set_immutable()
        return g


######################### decoders #################################

###### TamoBargIterariveDecoder #######
class TamoBargIterariveDecoder(Decoder):

    def __init__(self, code, max_num_of_itr=10):
        super().__init__(code, code.ambient_space(), "VectorEncoder")
        self._max_num_of_itr = max_num_of_itr
        C = self.code()
        self._F = C.base_field()
        self._partitions, _ = C.partition()
        self._evalpts_idx = C.evaluation_points_idx()

        i = 0
        self._partitions_local_decoders = []
        self._partitions_local_codes = []
        for partition in self._partitions:
            code_tmp = []
            decoder_tmp = []
            for coset in partition:
                local_grs = GeneralizedReedSolomonCode(coset, C.locality()[i])
                code_tmp.append(local_grs)
                decoder_tmp.append(local_grs.decoder("Gao"))

            self._partitions_local_decoders.append(decoder_tmp)
            self._partitions_local_codes.append(code_tmp)
            i += 1

    def _repr_(self):
        return "Iterative Decoder for %s" % self.code()

    def decode_to_code(self, r):
        uncorr_err = True
        num_itr = 0
        while ((uncorr_err is True) & (num_itr < self._max_num_of_itr)):
            uncorrectable_err = [False, False]
            for i in range(0, len(self._partitions)):
                partition = self._partitions[i]
                for j in range(0, len(partition)):
                    r_list = []
                    for k in range(0, len(partition[j])):
                        r_list.append(r[self._evalpts_idx[i][j][k]])

                    try:
                        corr_c = self._partitions_local_decoders[i][j].decode_to_code(
                            vector(self._F, r_list))
                    except:  # Decoding fails
                        uncorrectable_err[i] = True
                        corr_c = vector(self._F, r_list)

                    # Fix error in r
                    for l in range(0, len(partition[j])):
                        r[self._evalpts_idx[i][j][l]] = corr_c[l]

            if (uncorrectable_err[0] == False & uncorrectable_err[1] == False):
                uncorr_err = False
            num_itr += 1

        return r, num_itr

###### Error-and-Erasure Iterarive Decoder #######


class TamoBargIterariveEEDecoder(Decoder):

    def __init__(self, code, max_num_of_itr = 10, erasure_start_itr = None):
        super().__init__(code, code.ambient_space(), "VectorEncoder")
        self._max_num_of_itr = max_num_of_itr
        if erasure_start_itr is None:
            self._erasure_start_itr = self._max_num_of_itr // 2
        else:
            self._erasure_start_itr = erasure_start_itr

        C = self.code()
        self._F = C.base_field()
        self._partitions, _ = C.partition()
        self._evalpts_idx = C.evaluation_points_idx()

        i = 0
        self._partitions_local_decoders = []
        self._partitions_local_codes = []
        for partition in self._partitions:
            code_tmp = []
            decoder_tmp = []
            for coset in partition:
                local_grs = GeneralizedReedSolomonCode(coset, C.locality()[i])
                code_tmp.append(local_grs)
                decoder_tmp.append(local_grs.decoder("ErrorErasure"))

            self._partitions_local_decoders.append(decoder_tmp)
            self._partitions_local_codes.append(code_tmp)
            i += 1

    def _repr_(self):
        return "Iterative Decoder Error-and-Erasure for %s" % self.code()

    def decode_to_code(self, r):
        erasure_vector = zero_vector(GF(2), len(r))
        uncorr_err = True
        num_itr = 0
        while ((uncorr_err is True) & (num_itr < self._max_num_of_itr)):
            uncorrectable_err = [False, False]
            due_status = []
            for i in range(0, len(self._partitions)):
                partition = self._partitions[i]
                due_tmp = []
                for j in range(0, len(partition)):
                    coset = partition[j]
                    r_list = []
                    erasure_vector_list = []
                    for k in range(0, len(coset)):
                        r_list.append(r[self._evalpts_idx[i][j][k]])
                        erasure_vector_list.append(erasure_vector[self._evalpts_idx[i][j][k]])

                    word_and_erasure_vector = vector(self._F, r_list), vector(GF(2), erasure_vector_list)
                    erasure_vec = vector(GF(2), erasure_vector_list)

                    try:
                        corr_c = self._partitions_local_decoders[i][j].decode_to_code(word_and_erasure_vector)
                        due_tmp.append(False)
                    except:  # Decoding fails
                        due_tmp.append(True)
                        uncorrectable_err[i] = True
                        corr_c = vector(self._F, r_list)

                    # TODO: Implement single erasure decoder.
                    if erasure_vec.hamming_weight() == self._partitions_local_codes[i][j].minimum_distance()-1:
                        pts = []
                        err = []
                        for k in range(0, len(coset)):
                            if (erasure_vec[k] == 0):
                                pts.append((coset[k], r_list[k]))
                            else:
                                err.append(k)

                        R = PolynomialRing(self.code().base_field(), 'x')
                        d = R.lagrange_polynomial(pts)

                        corr_c = vector(self._F, r_list)
                        for e in err:
                            corr_c[e] = d(coset[e])

                    # Fix error in r
                    for l in range(0, len(coset)):
                        r[self._evalpts_idx[i][j][l]] = corr_c[l]

                due_status.append(due_tmp)

            if (uncorrectable_err[0] == False & uncorrectable_err[1] == False):
                uncorr_err = False

            num_itr += 1

            if (num_itr >= self._erasure_start_itr):
                # Loop over Partitions (Two partitions)
                for s in range(0, len(due_status)):
                    # Loop over Cosets
                    for t in range(0, len(self._partitions[s])):
                        # Loop over pts in coset
                        for m in range(0, len(self._partitions[s][t])):
                            # First partition loop marks erasure on all pts in coset with DUE.
                            if ((due_status[s][t] == True) & (s == 0)):
                                erasure_vector[self._evalpts_idx[s][t][m]] = 1
                            # Second parition loop cleans up erasure marking on pts belonging to coset with no errors.
                            elif (due_status[s][t] == False):
                                erasure_vector[self._evalpts_idx[s][t][m]] = 0

        corr_c_and_erasure_vector = r, vector(GF(2), erasure_vector)
        return corr_c_and_erasure_vector, num_itr

###### TamoBargGlobalDecoder #######


class TamoBargGlobalDecoder(Decoder):

    def __init__(self, code):
        super().__init__(code, code.ambient_space(), "VectorEncoder")
        self._global_grs = GeneralizedReedSolomonCode(
            code.evaluation_points(), code.length()-code.design_distance()+1)
        self._global_decoder = self._global_grs.decoder("Gao")

    def _repr_(self):
        return "Global Decoder for %s" % self.code()

    def decode_to_code(self, r):
        return self._global_decoder.decode_to_code(r)


###### TamoBargGlobal-Error-and-Erasure-Decoder #######
class TamoBargGlobalEEDecoder(Decoder):

    def __init__(self, code):
        super().__init__(code, code.ambient_space(), "VectorEncoder")
        self._global_grs = GeneralizedReedSolomonCode(
            code.evaluation_points(), code.length()-code.design_distance()+1)
        self._global_decoder = self._global_grs.decoder("ErrorErasure")

    def _repr_(self):
        return "Global Erasure-and-Error Decoder for %s" % self.code()

    def decode_to_code(self, word_and_erasure_vector):
        word, erasure_vector = word_and_erasure_vector

        # TODO: Implement single erasure decoder.
        if erasure_vector.hamming_weight() == self._global_grs.minimum_distance()-1:
            evalpts = self.code().evaluation_ponts()
            pts = []
            era = []
            for i in range(0, len(erasure_vector)):
                if (erasure_vector[i] == 0):
                    pts.append((evalpts[i], word[i]))
                else:  # Erasure idx
                    era.append(i)

            # Compute Encoding Poly
            R = PolynomialRing(self.code().base_field(), 'x')
            d = R.lagrange_polynomial(pts)

            corr_c = word
            # Compute value in erasure pts
            for e in era:
                corr_c[e] = d(evalpts[e])
        else:
            corr_c = self._global_decoder.decode_to_code(
                word_and_erasure_vector)

        return corr_c


###### TamoBargTwoStepDecoder #######
class TamoBargTwoStepDecoder(Decoder):

    def __init__(self, code, max_num_of_itr=10):
        super().__init__(code, code.ambient_space(), "VectorEncoder")
        self._iter_decoder = TamoBargIterariveDecoder(code, max_num_of_itr)
        self._global_decoder = TamoBargGlobalDecoder(code)

    def _repr_(self):
        return "Two-Steps Decoder for %s" % self.code()

    def decode_to_code(self, r):
        # First Step - Iterative Decoder
        tmp, num_itr = self._iter_decoder.decode_to_code(r)

        # Second Step - Global Decoder
        try:
            r = self._global_decoder.decode_to_code(tmp)
        except:
            r = tmp

        return r, num_itr


###### TamoBargTwoStepEEDecoder #######
class TamoBargTwoStepEEDecoder(Decoder):

    def __init__(self, code, max_num_of_itr=10, erasure_start_itr=None):
        super().__init__(code, code.ambient_space(), "VectorEncoder")
        self._iter_ee_decoder = TamoBargIterariveEEDecoder(code, max_num_of_itr, erasure_start_itr)
        self._global_ee_decoder = TamoBargGlobalEEDecoder(code)

    def _repr_(self):
        return "Two-Steps Erasure-and-Error Decoder for %s" % self.code()

    def decode_to_code(self, r):
        # First Step - Iterative Decoder
        tmp, num_itr = self._iter_ee_decoder.decode_to_code(r)
        num_of_rem_era = tmp[1].hamming_weight()

        # Second Step - Global Decoder
        try:
            r = self._global_ee_decoder.decode_to_code(tmp)
        except:
            r = tmp[0]

        return r, num_itr, num_of_rem_era


####################### registration ###############################
TamoBargCodeTwoSets._registered_encoders["VectorEncoder"] = TamoBargVectorEncoder
TamoBargCodeTwoSets._registered_decoders["IterativeDecoder"] = TamoBargIterariveDecoder
TamoBargCodeTwoSets._registered_decoders["TwoStepsDecoder"] = TamoBargTwoStepDecoder
TamoBargCodeTwoSets._registered_decoders["GlobalDecoder"] = TamoBargGlobalDecoder
TamoBargCodeTwoSets._registered_decoders["IterativeErasureErrorDecoder"] = TamoBargIterariveEEDecoder
TamoBargCodeTwoSets._registered_decoders["GlobalErasureErrorDecoder"] = TamoBargGlobalEEDecoder
TamoBargCodeTwoSets._registered_decoders["TwoStepsErasureErrorDecoder"] = TamoBargTwoStepEEDecoder
