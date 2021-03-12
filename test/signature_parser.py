from sig_utils import parse_signature, no_fwd_overload, no_rev_overload, ignored, eigen_types

class SignatureParser:
    def __init__(self, signature):
        self.return_type, self.function_name, self.stan_args = parse_signature(signature)
        self.name_counter = 0
        self.reset()

    def reset(self):
        self.code_list = []
    
    def number_arguments(self):
        return len(self.stan_args)

    def is_ode(self):
        return self.function_name.startswith("ode")
    
    def is_high_order(self):
        return any("=>" in arg for arg in self.stan_args)

    def is_rng(self):
        return self.function_name.endswith("_rng")
    
    def is_fwd_compatible(self):
        return not (self.function_name in no_fwd_overload or self.is_rng())

    def is_rev_compatible(self):
        return not (self.function_name in no_rev_overload or self.is_rng())

    def is_ignored(self):
        return self.function_name in ignored

    def has_vector_arg(self):
        return any(arg in eigen_types for arg in self.stan_args)

    def returns_int(self):
        return "int" in self.return_type