from signature_parser import SignatureParser
import pytest

@pytest.fixture
def add():
    return SignatureParser("add(real, vector) => vector")

@pytest.fixture
def ode():
    return SignatureParser("ode_adams((real, vector, ostream_ptr, vector) => vector, vector, real, array[] real, ostream_ptr, vector) => array[] vector")

@pytest.fixture
def rng():
    return SignatureParser("weibull_rng(row_vector, array[] real) => array[] real")

@pytest.fixture
def hmm():
    return SignatureParser("hmm_hidden_state_prob(matrix, matrix, vector) => matrix")

@pytest.fixture
def if_else():
    return SignatureParser("if_else(int, int, int) => int")

def test_number_arguments(add, ode):
    assert add.number_arguments() == 2
    assert ode.number_arguments() == 6

def test_ode(add, ode):
    assert add.is_ode() == False
    assert ode.is_ode() == True

def test_is_high_order(add, ode):
    assert add.is_high_order() == False
    assert ode.is_high_order() == True

def test_is_rng(add, rng):
    assert add.is_rng() == False
    assert rng.is_rng() == True

def test_is_fwd_compatible(add, ode):
    assert add.is_fwd_compatible() == True
    assert ode.is_fwd_compatible() == False

def test_is_rev_compatible(add, hmm):
    assert add.is_rev_compatible() == True
    assert hmm.is_rev_compatible() == False

def test_is_ignored(add, if_else):
    assert add.is_ignored() == False
    assert if_else.is_ignored() == True

def test_has_vector_arg(add, if_else):
    assert add.has_vector_arg() == True
    assert if_else.has_vector_arg() == False

def test_returns_int(add, if_else):
    assert add.returns_int() == False
    assert if_else.returns_int() == True