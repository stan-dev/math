from sig_utils import parse_signature, no_fwd_overload, no_rev_overload, ignored

class SignatureParser:
    """
    SignatureParser parses Stanc3 function signatures and returns helpful information about them

    :param signature: stanc3 function signature string
    """
    def __init__(self, signature):
        self.return_type, self.function_name, self.stan_args = parse_signature(signature)
    
    def number_arguments(self):
        """Get the number of arguments in this signature"""
        return len(self.stan_args)

    def is_ode(self):
        """Return true if this signature is an ODE function"""
        return self.function_name.startswith("ode")
    
    def is_high_order(self):
        """Return true if this signature is a higher order function"""
        return any("=>" in arg for arg in self.stan_args)

    def is_rng(self):
        """Return true if this signature is a random number generator"""
        return self.function_name.endswith("_rng")
    
    def is_fwd_compatible(self):
        """Return true if this signature should compatible with forward mode autodiff"""
        return not (self.function_name in no_fwd_overload or self.is_rng())

    def is_rev_compatible(self):
        """Return true if this signature should be compatible with reverse mode autodiff"""
        return not (self.function_name in no_rev_overload or self.is_rng())

    def is_ignored(self):
        """Return true if this signature should be ignored by any tests (they are functions in stanc3 not defined in Math)"""
        return self.function_name in ignored

    def has_eigen_compatible_arg(self):
        """Return true if any argument is vector-like (can be an Eigen c++ type)"""
        return any(arg in ("matrix", "vector", "row_vector") for arg in self.stan_args)

    def returns_int(self):
        """Return true if the function returns an int"""
        return "int" in self.return_type