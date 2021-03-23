from signature_parser import SignatureParser
import unittest

class SignatureParserTest(unittest.TestCase):
    def setUp(self):
        self.add = SignatureParser("add(real, vector) => vector")
        self.ode = SignatureParser("ode_adams((real, vector, ostream_ptr, vector) => vector, vector, real, array[] real, ostream_ptr, vector) => array[] vector")
        self.rng = SignatureParser("weibull_rng(row_vector, array[] real) => array[] real")
        self.hmm = SignatureParser("hmm_hidden_state_prob(matrix, matrix, vector) => matrix")
        self.if_else = SignatureParser("if_else(int, int, int) => int")

    def test_number_arguments(self):
        self.assertEquals(self.add.number_arguments(), 2)
        self.assertEquals(self.ode.number_arguments(), 6)

    def test_ode(self):
        self.assertEquals(self.add.is_ode(), False)
        self.assertEquals(self.ode.is_ode(), True)

    def test_is_high_order(self):
        self.assertEquals(self.add.is_high_order(), False)
        self.assertEquals(self.ode.is_high_order(), True)

    def test_is_rng(self):
        self.assertEquals(self.add.is_rng(), False)
        self.assertEquals(self.rng.is_rng(), True)

    def test_is_fwd_compatible(self):
        self.assertEquals(self.add.is_fwd_compatible(), True)
        self.assertEquals(self.ode.is_fwd_compatible(), False)

    def test_is_rev_compatible(self):
        self.assertEquals(self.add.is_rev_compatible(), True)
        self.assertEquals(self.hmm.is_rev_compatible(), False)

    def test_is_ignored(self):
        self.assertEquals(self.add.is_ignored(), False)
        self.assertEquals(self.if_else.is_ignored(), True)

    def test_has_vector_arg(self):
        self.assertEquals(self.add.has_eigen_compatible_arg(), True)
        self.assertEquals(self.if_else.has_eigen_compatible_arg(), False)

    def test_returns_int(self):
        self.assertEquals(self.add.returns_int(), False)
        self.assertEquals(self.if_else.returns_int(), True)

if __name__ == '__main__':
    unittest.main()