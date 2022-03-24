import collections
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
        """
        Add a statement to the code generator

        :param statement: An object of type statement_types.CppStatement
        """
        if not isinstance(statement, statement_types.CppStatement):
            raise TypeError(
                "Argument to FunctionGenerator._add_statement must be an instance of an object that inherits from CppStatement"
            )

        self.code_list.append(statement)
        return statement

    def _get_next_name_suffix(self):
        """Get the next available"""
        self.name_counter += 1
        return repr(self.name_counter - 1)

    def cpp(self):
        """Generate and return the c++ code corresponding to the list of statements in the code generator"""
        return os.linesep.join(statement.cpp() for statement in self.code_list)

    def build_arguments(self, signature_parser, arg_overloads, size):
        """
        Generate argument variables for each of the arguments in the given signature_parser
        with the given overloads in arg_overloads and with the given size

        :param signature_parser: An instance of SignatureParser
        :param arg_overloads: A list of argument overloads (Prim/Fwd/Rev/etc.) as strings
        :param size: Size of matrix-like arguments. This is not used for array arguments (which will effectively all be size 1)
        """
        arg_list = []
        for n, (overload, stan_arg) in enumerate(
            zip(arg_overloads, signature_parser.stan_args)
        ):
            suffix = self._get_next_name_suffix()

            number_nested_arrays, inner_type = parse_array(stan_arg)

            # Check if argument is differentiable
            if inner_type == "int" or n in non_differentiable_args.get(
                signature_parser.function_name, []
            ):
                overload = "Prim"

            # By default the variable value is None and a default will be substituted
            value = None

            # Check for special arguments (constrained variables or types)
            try:
                special_arg = special_arg_values[signature_parser.function_name][n]
                if isinstance(special_arg, str):
                    inner_type = special_arg
                elif special_arg is not None:
                    value = special_arg
            except KeyError:
                pass

            # The first case here is used for the array initializers in sig_utils.special_arg_values
            # Everything else uses the second case
            if number_nested_arrays > 0 and isinstance(value, collections.Iterable):
                arg = statement_types.ArrayVariable(
                    overload,
                    "array" + suffix,
                    number_nested_arrays,
                    inner_type,
                    size=1,
                    value=value,
                )
            else:
                if inner_type == "int":
                    arg = statement_types.IntVariable("int" + suffix, value)
                elif inner_type == "real":
                    arg = statement_types.RealVariable(overload, "real" + suffix, value)
                elif inner_type == "complex":
                    arg = statement_types.ComplexVariable(
                        overload, "complex" + suffix, value
                    )
                elif inner_type in (
                    "vector",
                    "row_vector",
                    "matrix",
                    "complex_vector",
                    "complex_row_vector",
                    "complex_matrix",
                ):
                    arg = statement_types.MatrixVariable(
                        overload, "matrix" + suffix, inner_type, size, value
                    )
                elif inner_type == "rng":
                    arg = statement_types.RngVariable("rng" + suffix)
                elif inner_type == "ostream_ptr":
                    arg = statement_types.OStreamVariable("ostream" + suffix)
                elif inner_type == "scalar_return_type":
                    arg = statement_types.ReturnTypeTVariable(
                        "ret_type" + suffix, *arg_list
                    )
                elif inner_type == "simplex":
                    arg = statement_types.SimplexVariable(
                        overload, "simplex" + suffix, size, value
                    )
                elif inner_type == "positive_definite_matrix":
                    arg = statement_types.PositiveDefiniteMatrixVariable(
                        overload, "positive_definite_matrix" + suffix, size, value
                    )
                elif (
                    inner_type
                    == "(vector, vector, data array[] real, data array[] int) => vector"
                ):
                    arg = statement_types.AlgebraSolverFunctorVariable(
                        "functor" + suffix
                    )
                elif inner_type == "(real, vector, ostream_ptr, vector) => vector":
                    arg = statement_types.OdeFunctorVariable("functor" + suffix)
                else:
                    raise Exception("Inner type " + inner_type + " not supported")

                if number_nested_arrays > 0:
                    self._add_statement(arg)
                    arg = statement_types.ArrayVariable(
                        overload,
                        "array" + suffix,
                        number_nested_arrays,
                        inner_type,
                        size=1,
                        value=arg,
                    )

            arg_list.append(self._add_statement(arg))

        if signature_parser.is_rng():
            arg_list.append(
                self._add_statement(
                    statement_types.RngVariable("rng" + self._get_next_name_suffix())
                )
            )

        return arg_list

    def add(self, arg1, arg2):
        """
        Generate code for arg1 + arg2

        :param arg1: First argument
        :param arg1: Second argument
        """
        return self._add_statement(
            statement_types.FunctionCall(
                "stan::math::add",
                "sum_of_sums" + self._get_next_name_suffix(),
                arg1,
                arg2,
            )
        )

    def convert_to_expression(self, arg, size=None):
        """
        Generate code to convert arg to an expression type of given size. If size is None, use the argument size

        :param arg: Argument to convert to expression
        """
        return self._add_statement(
            statement_types.ExpressionVariable(
                arg.name + "_expr" + self._get_next_name_suffix(), arg, size
            )
        )

    def expect_adj_eq(self, arg1, arg2):
        """
        Generate code that checks that the adjoints of arg1 and arg2 are equal

        :param arg1: First argument
        :param arg2: Second argument
        """
        return self._add_statement(
            statement_types.FunctionCall("stan::test::expect_adj_eq", None, arg1, arg2)
        )

    def expect_eq(self, arg1, arg2):
        """
        Generate code that checks that values of arg1 and arg2 are equal

        :param arg1: First argument
        :param arg2: Second argument
        """
        return self._add_statement(
            statement_types.FunctionCall("EXPECT_STAN_EQ", None, arg1, arg2)
        )

    def expect_leq_one(self, arg):
        """
        Generate code to check that arg is less than or equal to one

        :param arg: Argument to check
        """
        one = self._add_statement(
            statement_types.IntVariable("int" + self._get_next_name_suffix(), 1)
        )
        return self._add_statement(
            statement_types.FunctionCall("EXPECT_LE", None, arg, one)
        )

    def function_call_assign(self, cpp_function_name, *args):
        """
        Generate code to call the c++ function given by cpp_function_name with given args and assign the result to another variable

        :param cpp_function_name: c++ function name to call
        :param args: list of arguments to pass to function
        """
        return self._add_statement(
            statement_types.FunctionCall(
                cpp_function_name, "result" + self._get_next_name_suffix(), *args
            )
        )

    def grad(self, arg):
        """
        Generate code to call stan::test::grad(arg) (equivalent of arg.grad())

        :param arg: Argument to call grad on
        """
        return self._add_statement(
            statement_types.FunctionCall("stan::test::grad", None, arg)
        )

    def recover_memory(self):
        """Generate code to call stan::math::recover_memory()"""
        return self._add_statement(
            statement_types.FunctionCall("stan::math::recover_memory", None)
        )

    def recursive_sum(self, arg):
        """
        Generate code that repeatedly sums arg until all that is left is a scalar

        :param arg: Argument to sum
        """
        return self._add_statement(
            statement_types.FunctionCall(
                "stan::test::recursive_sum",
                "summed_result" + self._get_next_name_suffix(),
                arg,
            )
        )

    def to_var_value(self, arg):
        """
        Generate code to convert arg to a varmat

        :param arg: Argument to convert to varmat
        """
        return self._add_statement(
            statement_types.FunctionCall(
                "stan::math::to_var_value",
                arg.name + "_varmat" + self._get_next_name_suffix(),
                arg,
            )
        )
