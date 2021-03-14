import numbers
import os
import statement_types
from sig_utils import parse_array, non_differentiable_args, special_arg_values

class CodeGenerator:
    def __init__(self):
        self.name_counter = 0
        self.code_list = []

    def add_statement(self, statement):
        if not isinstance(statement, statement_types.CppStatement):
            raise Exception("Argument to FunctionGenerator.add_statement must be an instance of an object that inherits from CppStatement")

        self.code_list.append(statement)
        return statement

    def cpp(self):
        return os.linesep.join(statement.cpp() for statement in self.code_list)

    def get_next_name_suffix(self):
        self.name_counter += 1
        return repr(self.name_counter - 1)

    def build_arguments(self, signature_parser, arg_overloads, max_size = None):
        if not max_size:
            max_size = 2

        arg_list = []
        for n, (overload, stan_arg) in enumerate(zip(arg_overloads, signature_parser.stan_args)):
            suffix = self.get_next_name_suffix()

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

            if value is None or isinstance(value, numbers.Number):
                if inner_type == "int":
                    arg = statement_types.IntArgument("int" + suffix, value)
                elif inner_type == "real":
                    arg = statement_types.RealArgument(overload, "real" + suffix, value)
                elif inner_type in ("vector", "row_vector", "matrix"):
                    arg = statement_types.MatrixArgument(overload, "matrix" + suffix, inner_type, value)
                elif inner_type == "rng":
                    arg = statement_types.RngArgument("rng" + suffix)
                elif inner_type == "ostream_ptr":
                    arg = statement_types.OStreamArgument("ostream" + suffix)
                elif inner_type == "scalar_return_type":
                    arg = statement_types.ReturnTypeTArgument("ret_type" + suffix, *arg_list)
                elif inner_type == "simplex":
                    arg = statement_types.SimplexArgument(overload, "simplex" + suffix, value)
                elif inner_type == "positive_definite_matrix":
                    arg = statement_types.PositiveDefiniteMatrixArgument(overload, "positive_definite_matrix" + suffix, value)
                elif inner_type == "(vector, vector, data array[] real, data array[] int) => vector":
                    arg = statement_types.AlgebraSolverFunctorArgument("functor" + suffix)
                elif inner_type == "(real, vector, ostream_ptr, vector) => vector":
                    arg = statement_types.OdeFunctorArgument("functor" + suffix)
                else:
                    raise Exception("Inner type " + inner_type + " not supported")

                if number_nested_arrays > 0:
                    self.add_statement(arg)
                    arg = statement_types.ArrayArgument(overload, "array" + suffix, number_nested_arrays, inner_type, value = arg)
            else:
                arg = statement_types.ArrayArgument(overload, "array" + suffix, number_nested_arrays, inner_type, value = value)

            arg_list.append(self.add_statement(arg))

        if signature_parser.is_rng():
            arg_list.append(self.add_statement(statement_types.RngArgument("rng" + self.get_next_name_suffix())))

        return arg_list

    def add(self, arg1, arg2):
        return self.add_statement(statement_types.FunctionCallAssign("stan::math::add", "sum_of_sums" + self.get_next_name_suffix(), arg1, arg2))
    
    def convert_to_expression(self, arg):
        return self.add_statement(statement_types.ExpressionArgument(arg.overload, arg.name + "_expr" + self.get_next_name_suffix(), arg))
    
    def expect_adj_eq(self, arg1, arg2):
        return self.add_statement(statement_types.FunctionCall("stan::test::expect_adj_eq", arg1, arg2))

    def expect_eq(self, arg1, arg2):
        return self.add_statement(statement_types.FunctionCall("EXPECT_STAN_EQ", arg1, arg2))
    
    def expect_leq_one(self, arg):
        return self.add_statement(statement_types.FunctionCall("EXPECT_LEQ_ONE", arg))
    
    def function_call_assign(self, cpp_function_name, *args):
        return self.add_statement(statement_types.FunctionCallAssign(cpp_function_name, "result" + self.get_next_name_suffix(), *args))

    def grad(self, arg):
        return self.add_statement(statement_types.FunctionCall("stan::test::grad", arg))
    
    def recover_memory(self):
        return self.add_statement(statement_types.FunctionCall("stan::math::recover_memory"))

    def recursive_sum(self, arg):
        return self.add_statement(statement_types.FunctionCallAssign("stan::test::recursive_sum", "summed_result" + self.get_next_name_suffix(), arg))
    
    def to_var_value(self, arg):
        return self.add_statement(statement_types.FunctionCallAssign("stan::math::to_var_value", arg.name + "_varmat" + self.get_next_name_suffix(), arg))
