import numbers
import os
import statement_types
from sig_utils import parse_array, non_differentiable_args, special_arg_values

class CodeGenerator:
    """
    This class generates C++ to test Stan functions
    """
    def __init__(self):
        self.name_counter = 0
        self.code_list = []

    def _add_statement(self, statement):
        """Add a statement to the code generator"""
        if not isinstance(statement, statement_types.CppStatement):
            raise Exception("Argument to FunctionGenerator._add_statement must be an instance of an object that inherits from CppStatement")

        self.code_list.append(statement)
        return statement

    def _get_next_name_suffix(self):
        """Get the next available """
        self.name_counter += 1
        return repr(self.name_counter - 1)

    def cpp(self):
        """Generate and return the c++ code corresponding to the list of statements in the code generator"""
        return os.linesep.join(statement.cpp() for statement in self.code_list)

    def build_arguments(self, signature_parser, arg_overloads, size = None):
        """
        Generate argument variables for each of the arguments in the given signature_parser
        with the given overloads in arg_overloads and with the given size

        :param signature_parser: An instance of SignatureParser
        :param arg_overloads: A list of argument overloads (Prim/Fwd/Rev/etc.) as strings
        :param size: Size of matrix-like arguments. This is not used for array arguments (which will effectively all be size 1)
        """
        arg_list = []
        for n, (overload, stan_arg) in enumerate(zip(arg_overloads, signature_parser.stan_args)):
            suffix = self._get_next_name_suffix()

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
                    arg = statement_types.IntVariable("int" + suffix, value)
                elif inner_type == "real":
                    arg = statement_types.RealVariable(overload, "real" + suffix, value)
                elif inner_type in ("vector", "row_vector", "matrix"):
                    arg = statement_types.MatrixVariable(overload, "matrix" + suffix, inner_type, value, size)
                elif inner_type == "rng":
                    arg = statement_types.RngVariable("rng" + suffix)
                elif inner_type == "ostream_ptr":
                    arg = statement_types.OStreamVariable("ostream" + suffix)
                elif inner_type == "scalar_return_type":
                    arg = statement_types.ReturnTypeTVariable("ret_type" + suffix, *arg_list)
                elif inner_type == "simplex":
                    arg = statement_types.SimplexVariable(overload, "simplex" + suffix, value, size)
                elif inner_type == "positive_definite_matrix":
                    arg = statement_types.PositiveDefiniteMatrixVariable(overload, "positive_definite_matrix" + suffix, value, size)
                elif inner_type == "(vector, vector, data array[] real, data array[] int) => vector":
                    arg = statement_types.AlgebraSolverFunctorVariable("functor" + suffix)
                elif inner_type == "(real, vector, ostream_ptr, vector) => vector":
                    arg = statement_types.OdeFunctorVariable("functor" + suffix)
                else:
                    raise Exception("Inner type " + inner_type + " not supported")

                if number_nested_arrays > 0:
                    self._add_statement(arg)
                    arg = statement_types.ArrayVariable(overload, "array" + suffix, number_nested_arrays, inner_type, arg, 1)
            else:
                arg = statement_types.ArrayVariable(overload, "array" + suffix, number_nested_arrays, inner_type, value, 1)

            arg_list.append(self._add_statement(arg))

        if signature_parser.is_rng():
            arg_list.append(self._add_statement(statement_types.RngVariable("rng" + self._get_next_name_suffix())))

        return arg_list

    def add(self, arg1, arg2):
        """Generate code for arg1 + arg2"""
        return self._add_statement(statement_types.FunctionCallAssign("stan::math::add", "sum_of_sums" + self._get_next_name_suffix(), arg1, arg2))
    
    def convert_to_expression(self, arg, size = None):
        """Generate code to convert arg to an expression type of given size. If size is None, use the argument size"""
        if not size:
            size = arg.size

        return self._add_statement(statement_types.ExpressionVariable(arg.name + "_expr" + self._get_next_name_suffix(), arg, size))
    
    def expect_adj_eq(self, arg1, arg2):
        """Generate code that checks that the adjoints of arg1 and arg2 are equal"""
        return self._add_statement(statement_types.FunctionCall("stan::test::expect_adj_eq", arg1, arg2))

    def expect_eq(self, arg1, arg2):
        """Generate code that checks that values of arg1 and arg2 are equal"""
        return self._add_statement(statement_types.FunctionCall("EXPECT_STAN_EQ", arg1, arg2))
    
    def expect_leq_one(self, arg):
        """Generate code to check that arg is less than or equal to one"""
        return self._add_statement(statement_types.FunctionCall("EXPECT_LEQ_ONE", arg))
    
    def function_call_assign(self, cpp_function_name, *args):
        """Generate code to call the c++ function given by cpp_function_name with given args and assign the result to another variable"""
        return self._add_statement(statement_types.FunctionCallAssign(cpp_function_name, "result" + self._get_next_name_suffix(), *args))

    def grad(self, arg):
        """Generate code to call stan::math::grad()"""
        return self._add_statement(statement_types.FunctionCall("stan::test::grad", arg))
    
    def recover_memory(self):
        """Generate code to call stan::math::recover_memory()"""
        return self._add_statement(statement_types.FunctionCall("stan::math::recover_memory"))

    def recursive_sum(self, arg):
        """Generate code that repeatedly sums arg until all that is left is a scalar"""
        return self._add_statement(statement_types.FunctionCallAssign("stan::test::recursive_sum", "summed_result" + self._get_next_name_suffix(), arg))
    
    def to_var_value(self, arg):
        """Generate code to convert arg to a varmat"""
        return self._add_statement(statement_types.FunctionCallAssign("stan::math::to_var_value", arg.name + "_varmat" + self._get_next_name_suffix(), arg))
