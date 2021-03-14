import numbers
import os

from sig_utils import overload_scalar, get_cpp_type, arg_types

class CppStatement:
    def __init__(self):
        raise Exception("CppStatement should never be instantiated directly")

    def is_reverse_mode(self):
        return False
    
    def is_matrix_like(self):
        return False

    def is_varmat_compatible(self):
        return False

    def is_expression(self):
        return False

class IntArgument(CppStatement):
    def __init__(self, name, value = None):
        self.name = name
        if value is None:
            self.value = 1
        else:
            self.value = value

    def cpp(self):
        return "int {name} = {value};".format(name = self.name, value = self.value)

class RealArgument(CppStatement):
    def __init__(self, overload, name, value = None):
        self.overload = overload
        self.name = name
        if value is None:
            self.value = 0.4
        else:
            self.value = value
    
    def is_reverse_mode(self):
        return self.overload.startswith("Rev")
    
    def cpp(self):
        scalar = overload_scalar[self.overload]

        return "{type} {name} = {value};".format(type = scalar, name = self.name, value = self.value)

class MatrixArgument(CppStatement):
    def __init__(self, overload, name, stan_arg, value = None):
        self.overload = overload
        self.name = name
        self.stan_arg = stan_arg
        if value is None:
            self.value = 0.4
        else:
            self.value = value
    
    def is_reverse_mode(self):
        return self.overload.startswith("Rev")

    def is_matrix_like(self):
        return True
    
    def is_varmat_compatible(self):
        return True
    
    def is_varmat(self):
        return self.overload.endswith("Varmat")
    
    def cpp(self):
        scalar = overload_scalar[self.overload]
        cpp_arg_template = get_cpp_type(self.stan_arg)
        arg_type = cpp_arg_template.replace("SCALAR", scalar)

        return "auto {name} = stan::test::make_arg<{type}>({value}, 1);".format(name = self.name, type = arg_type, value = self.value)

class SimplexArgument(CppStatement):
    def __init__(self, overload, name, value = None):
        self.overload = overload
        self.name = name
        self.stan_arg = "vector"
        if value is None:
            self.value = 0.4
        else:
            self.value = value
    
    def is_reverse_mode(self):
        return self.overload.startswith("Rev")

    def is_matrix_like(self):
        return True
    
    def is_varmat_compatible(self):
        return True
    
    def cpp(self):
        scalar = overload_scalar[self.overload]
        cpp_arg_template = get_cpp_type(self.stan_arg)
        arg_type = cpp_arg_template.replace("SCALAR", scalar)

        return "auto {name} = stan::test::make_simplex<{type}>({value}, 1);".format(name = self.name, type = arg_type, value = self.value)

class PositiveDefiniteMatrixArgument(CppStatement):
    def __init__(self, overload, name, value = None):
        self.overload = overload
        self.name = name
        self.stan_arg = "matrix"
        if value is None:
            self.value = 0.4
        else:
            self.value = value
    
    def is_reverse_mode(self):
        return self.overload.startswith("Rev")

    def is_matrix_like(self):
        return True
    
    def is_varmat_compatible(self):
        return True
    
    def cpp(self):
        scalar = overload_scalar[self.overload]
        cpp_arg_template = get_cpp_type(self.stan_arg)
        arg_type = cpp_arg_template.replace("SCALAR", scalar)

        return "auto {name} = stan::test::make_pos_definite_matrix<{type}>({value}, 1);".format(name = self.name, type = arg_type, value = self.value)

# These are used with algebra_solver
class AlgebraSolverFunctorArgument(CppStatement):
    def __init__(self, name):
        self.name = name
    
    def cpp(self):
        return "stan::test::simple_eq_functor {name};".format(name = self.name)

# These are used with ode_adams
class OdeFunctorArgument(CppStatement):
    def __init__(self, name):
        self.name = name
    
    def cpp(self):
        return "stan::test::test_functor {name};".format(name = self.name)

class ReturnTypeTArgument(CppStatement):
    def __init__(self, name, *args):
        self.name = name
        arg_names = ["decltype({name})".format(name = arg.name) for arg in args]
        self.type_str = "stan::return_type_t<{names}>".format(names = ','.join(arg_names))
        self.__is_reverse_mode = any(arg.is_reverse_mode() for arg in args)
    
    def is_reverse_mode(self):
        return self.__is_reverse_mode
    
    def cpp(self):
        return "{type} {name} = 0;".format(type = self.type_str, name = self.name)

class RngArgument(CppStatement):
    def __init__(self, name):
        self.name = name
    
    def cpp(self):
        return "std::minstd_rand {name};".format(name = self.name)

class OStreamArgument(CppStatement):
    def __init__(self, name):
        self.name = name
    
    def cpp(self):
        return "std::ostream* {name} = &std::cout;".format(name = self.name)

class ArrayArgument(CppStatement):
    def __init__(self, overload, name, number_nested_arrays, inner_type, value = None, size = None):
        self.overload = overload
        self.name = name
        self.number_nested_arrays = number_nested_arrays
        self.inner_type = inner_type
        self.value = value

        if size:
            self.size = size
        else:
            self.size = 1
    
    def is_reverse_mode(self):
        try:
            return self.value.is_reverse_mode()
        except AttributeError:
            return self.overload.startswith("Rev") and self.inner_type != "int"

    def is_varmat_compatible(self):
        try:
            return self.value.is_varmat_compatible()
        except AttributeError:
            return False
    
    def cpp(self):
        try:
            lhs_string = "decltype({name})".format(name = self.value.name)
            rhs_string = self.value.name
            for n in range(self.number_nested_arrays):
                lhs_string = "<" + lhs_string + ">"
                rhs_string = "{" + ",".join([rhs_string] * self.size) + "}"
            arg_type = "std::vector" + lhs_string

            return "{type} {name} = {rhs_string};".format(type = arg_type, name = self.name, lhs_string = lhs_string, rhs_string = rhs_string)
        except AttributeError:
            lhs_string = arg_types[self.inner_type]
            for n in range(self.number_nested_arrays):
                lhs_string = "<" + lhs_string + ">"
            arg_type = "std::vector" + lhs_string

            value_string = repr(self.value).replace('[', '{').replace(']', '}').replace('(', '{').replace(')', '}')
            return "{type} {name} = {value};".format(type = arg_type, name = self.name, value = value_string)

class FunctionCallAssign(CppStatement):
    def __init__(self, function_name, name, *args):
        self.function_name = function_name
        self.name = name
        self.arg_str = ','.join(arg.name for arg in args)
        # This isn't always accurate cause the output could be an int or a boolean
        # or something even if the input is reverse mode
        self.__is_reverse_mode = any(arg.is_reverse_mode() for arg in args)

    def is_reverse_mode(self):
        return self.__is_reverse_mode

    def is_matrix_like(self):
        raise Exception("is_matrix_like not implemented for FunctionCallAssign")

    def is_varmat(self):
        raise Exception("is_varmat not implemented for FunctionCallAssign")

    def cpp(self):
        if self.name:
            return "auto {name} = stan::math::eval({function_name}({arg_str}));".format(name = self.name, function_name = self.function_name, arg_str = self.arg_str)
        else:
            return "{function_name}({arg_str});".format(function_name = self.function_name, arg_str = self.arg_str)

class FunctionCall(FunctionCallAssign):
    def __init__(self, function_name, *args):
        FunctionCallAssign.__init__(self, function_name, None, *args)

class ExpressionArgument(CppStatement):
    def __init__(self, overload, name, arg):
        self.overload = overload
        self.name = name
        self.counter = IntArgument("{name}_counter".format(name = self.name), 0)
        self.arg = arg
    
    def is_reverse_mode(self):
        return self.arg.is_reverse_mode()

    def is_matrix_like(self):
        return True
    
    def is_expression(self):
        return True

    def cpp(self):
        scalar = overload_scalar[self.overload]
        counter_op_name = "{name}_counter_op".format(name = self.name)

        code = (
            self.counter.cpp() + os.linesep +
            "stan::test::counterOp<{scalar}> {counter_op_name}(&{counter});".format(scalar = scalar, counter_op_name = counter_op_name, counter = self.counter.name) + os.linesep
        )

        if self.arg.stan_arg == "matrix":
            return code + "auto {name} = {arg}.block(0,0,1,1).unaryExpr({counter_op});".format(name = self.name, arg = self.arg.name, counter_op = counter_op_name)
        elif self.arg.stan_arg in ("vector", "row_vector"):
            return code + "auto {name} = {arg}.segment(0,1).unaryExpr({counter_op});".format(name = self.name, arg = self.arg.name, counter_op = counter_op_name)
        else:
            raise Exception("Can't make an expression out of a " + self.arg.stan_arg)
