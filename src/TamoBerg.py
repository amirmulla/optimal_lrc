from sage.all import *
from sage.coding.linear_code import AbstractLinearCode
from sage.coding.encoder import Encoder
from sage.coding.decoder import Decoder, DecodingError
from sage.coding.grs_code import GeneralizedReedSolomonCode

####################### code ###############################


class TamoBergCode(AbstractLinearCode):
    _registered_encoders = {}
    _registered_decoders = {}

    def __init__(self, base_field, length, dimension, locality, local_minimum_distance=2, sub_group_type="any"):

        if not base_field.is_finite():
            raise ValueError("base_field has to be a finite field")

        super().__init__(base_field, length, "VectorEncoder", "ErasureDecoder")

        self._dimension = dimension
        self._locality = locality
        self._local_minimum_distance = local_minimum_distance
        self._base_mult_group = base_field.list()[1:base_field.cardinality()]
        self._mult_sub_groups = self._list_mult_subgroups()
        self._add_sub_groups = self._list_add_subgroups()

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
                    print("Can't find any subgroup of size: ",locality + local_minimum_distance - 1)

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

        self._parition_size = int(length /(locality + local_minimum_distance - 1))
        self._partition = self._list_sub_group_cosets()[:self._parition_size]
        self._evaluation_points = flatten(self._partition)

    def __eq__(self, other):
        return isinstance(other, TamoBergCode) \
            and self.length() == other.length() \
            and self.dimension() == other.dimension() \
            and self.locality() == other.locality()

    def _repr_(self):
        return "[%s, %s, %s] Tamo-Berg Code over GF(%s)"\
            % (self.length(), self.dimension(), self.locality(), self.base_field().cardinality())

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
        p = self.base_field().characteristic()
        for i in range(1, self.base_field().degree()):
            tmp = []
            F = GF(p**i, repr="int")
            for elm in F:
                tmp.append(F.fetch_int(int(str(elm))))

            tmp.sort()
            additive_subgroups[len(tmp)] = tmp

        return additive_subgroups

    def _list_sub_group_cosets(self):
        F = self.base_field()
        h_cosets = []
        for g in self._base_mult_group:
            tmp = []
            for h in self._sub_group:
                if (self._sub_group_type == "mult"):
                    tmp.append(h * g)
                elif (self._sub_group_type == "add"):
                    tmp.append(h + g)
            tmp.sort()
            h_cosets.append(tmp)

        return list(map(list, set(map(lambda i: tuple(i), h_cosets))))

    def partition(self):
        return self._partition, self._parition_size

    def sub_group(self):
        return self._sub_group, self._sub_group_type

    def evaluation_points(self):
        return self._evaluation_points

    def locality(self):
        return self._locality

    def global_minimum_distance(self):
        return (self.length() - self.dimension() -
                ceil(self.dimension() / self.locality()) + 2)

    def local_minimum_distance(self):
        return self._local_minimum_distance


####################### encoders ###############################


class TamoBergVectorEncoder(Encoder):

    def __init__(self, code):
        self._partition, self._partition_size = code.partition()
        self._sub_group, self._sub_group_type = code.sub_group()
        self._base_poly_ring = PolynomialRing(code.base_field(), 'x')
        self._annihilator = self._calc_annihilator(code.evaluation_points(),self._base_poly_ring)
        self._quotient_poly_ring = self._base_poly_ring.quotient(self._annihilator, 'x')
        self._good_poly = self._calc_good_poly(self._sub_group, self._sub_group_type, self._quotient_poly_ring)
        self._algebra_basis = self._calc_algebra_basis(self._good_poly, self._partition_size)
        self._enc_basis = self._calc_enc_basis(self._algebra_basis, code.locality(), self._quotient_poly_ring)

        self._max_dimension = len(self._enc_basis)
        if(code.dimension() > self._max_dimension):
            raise ValueError("code dimension is too big. maximum code dimension: ", self._max_dimension)
        
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
            good_poly = 1
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

    def enc_basis(self):
        return self._enc_basis

    def enc_poly(self):
        return self._enc_poly

    @cached_method
    def generator_matrix(self):
        C = self.code()
        alphas = C.evaluation_points()
        g = matrix(C.base_field(), C.dimension(), C.length(),
                   lambda i, j: self._enc_basis[i].lift()(alphas[j]))
        g.set_immutable()
        return g


######################### decoders #################################


class TamoBergErasureDecoder(Decoder):

    def __init__(self, code):
        input_space = cartesian_product([code.ambient_space(),VectorSpace(GF(2),code.ambient_space().dimension())])
        super().__init__(code, input_space, "VectorEncoder")

    def _repr_(self):
        return "Erasure Decoder for %s" % self.code()

    def decode_to_code(self, word_and_erasure_vector):
        C = self.code()
        word, erasure_vector = word_and_erasure_vector
        n, k = C.length(), C.dimension()
        m = C.locality() + C.local_minimum_distance() - 1
        erasure_vector_list = [erasure_vector[x:x + m] for x in range(0, len(erasure_vector), m)]
        word_list = [word[x:x + m] for x in range(0, len(word), m)]
        partition, _ = C.partition()

        for i in range(0, len(erasure_vector_list)):
            erasure_vec = erasure_vector_list[i]
            if erasure_vec.hamming_weight() >= self.code().local_minimum_distance():
                raise DecodingError("Too many erasures in the received word")
            elif erasure_vec.hamming_weight() != 0:
                pts = []
                err = []
                for j in range(0, m):
                    if (erasure_vec[j] == 0):
                        pts.append((partition[i][j], word_list[i][j]))
                    else:
                        err.append(j)

                R = PolynomialRing(C.base_field(), 'x')
                d = R.lagrange_polynomial(pts)

                for e in err:
                    word_list[i][e] = d(partition[i][e])

        l = []
        for vec in word_list:
            for v in vec:
                l.append(v)

        return (vector(C.base_field(), l))


class TamoBergErrorDecoder(Decoder):

    def __init__(self, code):
        input_space = cartesian_product([code.ambient_space(),VectorSpace(GF(2),code.ambient_space().dimension())])
        super().__init__(code, input_space, "VectorEncoder")

    def _repr_(self):
        return "Erasure Decoder for %s" % self.code()

    def decode_to_code(self, r):
        C = self.code()
        F = C.base_field()
        n, k = C.length(), C.dimension()

        m = C.locality() + C.local_minimum_distance() - 1
        r_list = [r[x:x + m] for x in range(0, len(r), m)]

        partition, _ = C.partition()

        corr_c = []

        for i in range(0, len(r_list)):
            try:
                local_grs = GeneralizedReedSolomonCode(evaluation_points=partition[i], dimension=C.locality())
                local_decoder = local_grs.decoder("BerlekampWelch")
                local_corr_c = local_decoder.decode_to_code(vector(F, r_list[i]))
                corr_c.append(local_corr_c.list())
            except:
                corr_c.append(r_list[i].list())

        corr_c = flatten(corr_c)

        return vector(F, corr_c)


####################### registration ###############################
TamoBergCode._registered_encoders["VectorEncoder"] = TamoBergVectorEncoder
TamoBergCode._registered_decoders["ErasureDecoder"] = TamoBergErasureDecoder
TamoBergCode._registered_decoders["ErrorDecoder"] = TamoBergErrorDecoder
