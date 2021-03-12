import numbers
import os
from sig_utils import overload_scalar, get_cpp_type, arg_types, parse_array, non_differentiable_args, special_arg_values

class CppStatement:
    def __init__(self):
        raise Exception("CppStatement should never be instantiated directly")

    def is_reverse_mode(self):
        return False
    
    def is_matrix_like(self):
        return False
    
    def is_varmat(self):
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
        return f"int {self.name} = {self.value};"

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

        return f"{scalar} {self.name} = {self.value};"

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
    
    def is_varmat(self):
        return self.overload.endswith("Varmat")
    
    def cpp(self):
        scalar = overload_scalar[self.overload]
        cpp_arg_template = get_cpp_type(self.stan_arg)
        arg_type = cpp_arg_template.replace("SCALAR", scalar)

        return f"auto {self.name} = stan::test::make_arg<{arg_type}>({self.value}, 1);"

class SimplexArgument(CppStatement):
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
    
    def is_varmat(self):
        return self.overload.endswith("Varmat")
    
    def cpp(self):
        scalar = overload_scalar[self.overload]
        cpp_arg_template = get_cpp_type(self.stan_arg)
        arg_type = cpp_arg_template.replace("SCALAR", scalar)

        return f"auto {self.name} = stan::test::make_simplex<{arg_type}>({self.value}, 1);"

class PositiveDefiniteMatrixArgument(CppStatement):
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
    
    def is_varmat(self):
        return self.overload.endswith("Varmat")
    
    def cpp(self):
        scalar = overload_scalar[self.overload]
        cpp_arg_template = get_cpp_type(self.stan_arg)
        arg_type = cpp_arg_template.replace("SCALAR", scalar)

        return f"auto {self.name} = stan::test::make_pos_definite_matrix<{arg_type}>({self.value}, 1);"

# These are used with algebra_solver
class AlgebraSolverFunctorArgument(CppStatement):
    def __init__(self, name):
        self.name = name
    
    def cpp(self):
        return f"stan::test::simple_eq_functor {self.name};"

# These are used with ode_adams
class OdeFunctorArgument(CppStatement):
    def __init__(self, name):
        self.name = name
    
    def cpp(self):
        return f"stan::test::test_functor {self.name};"

class ReturnTypeTArgument(CppStatement):
    def __init__(self, name, *args):
        self.name = name
        arg_names = [f'decltype({arg.name})' for arg in args]
        self.type_str = f"stan::return_type_t<{','.join(arg_names)}>"
        self.__is_reverse_mode = any(arg.is_reverse_mode() for arg in args)
    
    def is_reverse_mode(self):
        return self.__is_reverse_mode
    
    def cpp(self):
        return f"{self.type_str} {self.name} = 0;"

class RngArgument(CppStatement):
    def __init__(self, name):
        self.name = name
    
    def cpp(self):
        return f"std::minstd_rand {self.name};"

class OStreamArgument(CppStatement):
    def __init__(self, name):
        self.name = name
    
    def cpp(self):
        return f"std::ostream* {self.name} = &std::cout;"

class ArrayArgument(CppStatement):
    def __init__(self, overload, name, number_nested_arrays, inner_type, value = None, size = None):
        self.overload = overload
        self.name = name
        self.number_nested_arrays = number_nested_arrays
        self.inner_type = inner_type

        if not value:
            if inner_type == "int":
                self.value = 1
            else:
                self.value = 0.4
        else:
            self.value = value

        if size:
            self.size = size
        else:
            self.size = 1
    
    def is_reverse_mode(self):
        return self.overload.startswith("Rev")
    
    def cpp(self):
        scalar = overload_scalar[self.overload]
        arg_type = arg_types[self.inner_type]
        for n in range(self.number_nested_arrays):
            arg_type = "std::vector<{}>".format(arg_type)
        arg_type = arg_type.replace("SCALAR", scalar)

        if isinstance(self.value, numbers.Number):
            value_string = f"{self.value}"
        else:
            value_string = f"{{{','.join([repr(v) for v in self.value])}}}"

        if isinstance(self.value, numbers.Number):
            return f"auto {self.name} = stan::test::make_arg<{arg_type}>({self.value}, {self.size});"
        else:
            value_string = repr(self.value).replace('[', '{').replace(']', '}').replace('(', '{').replace(')', '}')
            return f"{arg_type} {self.name} = {value_string};"

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
            return f"auto {self.name} = stan::math::eval({self.function_name}({self.arg_str}));"
        else:
            return f"{self.function_name}({self.arg_str});"

class FunctionCall(FunctionCallAssign):
    def __init__(self, function_name, *args):
        FunctionCallAssign.__init__(self, function_name, None, *args)

class ExpressionArgument(CppStatement):
    def __init__(self, overload, name, arg):
        self.overload = overload
        self.name = name
        self.counter = IntArgument(f"{self.name}_counter", 0)
        self.arg = arg
    
    def is_reverse_mode(self):
        return self.arg.is_reverse_mode()

    def is_matrix_like(self):
        return True
    
    def is_expression(self):
        return True

    def cpp(self):
        scalar = overload_scalar[self.overload]
        counter_op_name = f"{self.name}_counter_op"

        code = (
            self.counter.cpp() + os.linesep +
            f"stan::test::counterOp<{scalar}> {counter_op_name}(&{self.counter.name});" + os.linesep
        )

        if self.arg.stan_arg == "matrix":
            return code + f"auto {self.name} = {self.arg.name}.block(0,0,1,1).unaryExpr({counter_op_name});"
        elif self.arg.stan_arg in ("vector", "row_vector"):
            return code + f"auto {self.name} = {self.arg.name}.segment(0,1).unaryExpr({counter_op_name});"
        else:
            raise Exception(f"Can't make an expression out of a {self.arg.stan_arg}")

class CodeGenerator:
    def __init__(self):
        self.name_counter = 0
        self.code_list = []

    def add_op(self, statement):
        if not isinstance(statement, CppStatement):
            raise Exception("Argument to FunctionGenerator.add_op must be an instance of an object that inherits from CppStatement")

        self.code_list.append(statement)
        return statement

    def cpp(self):
        return os.linesep.join(statement.cpp() for statement in self.code_list)

    def get_next_name(self):
        self.name_counter += 1
        return repr(self.name_counter - 1)

    def build_arguments(self, signature_parser, arg_overloads, max_size = None):
        if not max_size:
            max_size = 2

        arg_list = []
        for overload, stan_arg in zip(arg_overloads, signature_parser.stan_args):
            n = self.get_next_name()

            number_nested_arrays, inner_type = parse_array(stan_arg)

            value = None

            # Check if argument is differentiable
            try:
                if inner_type == "int":
                    overload = "Prim"
                if n in non_differentiable_args[signature_parser.function_name]:
                    overload = "Prim"
            except KeyError:
                pass

            # Check for special arguments (constrained variables or types)
            try:
                special_arg = special_arg_values[signature_parser.function_name][n]
                if isinstance(special_arg, str):
                    inner_type = special_arg
                elif special_arg is not None:
                    value = special_arg
            except KeyError:
                pass

            if number_nested_arrays == 0:
                if inner_type == "int":
                    arg = IntArgument(f"int{n}", value)
                elif inner_type == "real":
                    arg = RealArgument(overload, f"real{n}", value)
                elif inner_type in ("vector", "row_vector", "matrix"):
                    arg = MatrixArgument(overload, f"matrix{n}", stan_arg, value)
                elif inner_type == "rng":
                    arg = RngArgument(f"rng{n}")
                elif inner_type == "ostream_ptr":
                    arg = OStreamArgument(f"ostream{n}")
                elif inner_type == "scalar_return_type":
                    arg = ReturnTypeTArgument(f"ret_type{n}", *arg_list)
                elif inner_type == "simplex":
                    arg = SimplexArgument(overload, f"simplex{n}", stan_arg, value)
                elif inner_type == "positive_definite_matrix":
                    arg = PositiveDefiniteMatrixArgument(overload, f"positive_definite_matrix{n}", stan_arg, value)
                elif inner_type == "(vector, vector, data array[] real, data array[] int) => vector":
                    arg = AlgebraSolverFunctorArgument(f"functor{n}")
                elif inner_type == "(real, vector, ostream_ptr, vector) => vector":
                    arg = OdeFunctorArgument(f"functor{n}")
                else:
                    raise Exception(f"Inner type '{inner_type}' not supported")
            else:
                arg = ArrayArgument(overload, f"array{n}", number_nested_arrays, inner_type, value = value)

            arg_list.append(self.add_op(arg))

        if signature_parser.is_rng():
            arg_list.append(self.add_op(RngArgument(f"rng{self.get_next_name()}")))

        return arg_list

    def add(self, arg1, arg2):
        return self.add_op(FunctionCallAssign("stan::math::add", f"sum_of_sums{self.get_next_name()}", arg1, arg2))
    
    def convert_to_expression(self, arg):
        return self.add_op(ExpressionArgument(arg.overload, f"{arg.name}_expr{self.get_next_name()}", arg))
    
    def expect_adj_eq(self, arg1, arg2):
        return self.add_op(FunctionCall("stan::test::expect_adj_eq", arg1, arg2))

    def expect_eq(self, arg1, arg2):
        return self.add_op(FunctionCall("EXPECT_STAN_EQ", arg1, arg2))
    
    def expect_leq_one(self, arg):
        return self.add_op(FunctionCall("EXPECT_LEQ_ONE", arg))
    
    def function_call_assign(self, cpp_function_name, *args):
        return self.add_op(FunctionCallAssign(cpp_function_name, f"result{self.get_next_name()}", *args))

    def grad(self, arg):
        return self.add_op(FunctionCall("stan::test::grad", arg))
    
    def recover_memory(self):
        return self.add_op(FunctionCall("stan::math::recover_memory"))

    def recursive_sum(self, arg):
        return self.add_op(FunctionCallAssign("stan::test::recursive_sum", f"summed_result{self.get_next_name()}", arg))
    
    def to_var_value(self, arg):
        return self.add_op(FunctionCallAssign("stan::math::to_var_value", f"{arg.name}_varmat{self.get_next_name()}", arg))
