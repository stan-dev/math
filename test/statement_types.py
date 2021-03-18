import numbers
import os

from sig_utils import overload_scalar, get_cpp_type, arg_types

class CppStatement:
    """Base class for all statements. Implements some default checks. Should not be used directly"""
    def __init__(self):
        raise Exception("CppStatement should never be instantiated directly")

    def is_reverse_mode(self):
        """By default, a statement is not reverse mode"""
        return False
    
    def is_eigen_compatible(self):
        """By default, a statement is not matrix-like"""
        return False

    def is_varmat_compatible(self):
        """By default, a statement is not compatible with varmat types"""
        return False

    def is_expression(self):
        """By default, a statement is not an expression"""
        return False

class IntVariable(CppStatement):
    """Represents an integer variable"""
    def __init__(self, name, value = None):
        """
        Initialize matrix
        
        :param overload: Type of overload as string (Prim/Fwd/Rev/etc.)
        :param name: C++ name of variable
        :param value: Scalar value to initialize matrix with
        """
        self.name = name
        if value is None:
            self.value = 1
        else:
            self.value = int(value)

    def cpp(self):
        """Generate C++"""
        return "int {name} = {value};".format(name = self.name, value = self.value)

class RealVariable(CppStatement):
    """Represents a scalar, real variable"""
    def __init__(self, overload, name, value = None):
        """
        Initialize matrix
        
        :param overload: Type of overload as string (Prim/Fwd/Rev/etc.)
        :param name: C++ name of variable
        :param value: Scalar value to initialize matrix with
        """
        self.overload = overload
        self.name = name
        if value is None:
            self.value = 0.4
        else:
            self.value = value
    
    def is_reverse_mode(self):
        """Return true if the overload is reverse mode"""
        return self.overload.startswith("Rev")
    
    def cpp(self):
        """Generate c++"""
        scalar = overload_scalar[self.overload]

        return "{type} {name} = {value};".format(type = scalar, name = self.name, value = self.value)

class MatrixVariable(CppStatement):
    """Represents a matrix variable"""
    def __init__(self, overload, name, stan_arg, size, value = None):
        """
        Initialize matrix
        
        :param overload: Type of overload as string (Prim/Fwd/Rev/etc.)
        :param name: C++ name of variable
        :param stan_arg: Stanc3 type string
        :param size: Number of rows/columns in matrix
        :param value: Scalar value to initialize matrix with
        """
        self.overload = overload
        self.name = name
        self.stan_arg = stan_arg
        self.size = size
        
        if value is None:
            self.value = 0.4
        else:
            self.value = value
    
    def is_reverse_mode(self):
        """Return true if the overload is reverse mode"""
        return self.overload.startswith("Rev")

    def is_eigen_compatible(self):
        """Return true (this is a matrix like variable)"""
        return True
    
    def is_varmat_compatible(self):
        """Return true (matrix-like types are varmat compatible)"""
        return True
    
    def cpp(self):
        """Generate C++"""
        scalar = overload_scalar[self.overload]
        cpp_arg_template = get_cpp_type(self.stan_arg)
        arg_type = cpp_arg_template.replace("SCALAR", scalar)

        return "auto {name} = stan::test::make_arg<{type}>({value}, {size});".format(name = self.name, type = arg_type, value = self.value, size = self.size)

class SimplexVariable(CppStatement):
    """Represents a simplex variable"""
    def __init__(self, overload, name, size, value = None):
        """
        Initialize simplex
        
        :param overload: Type of overload as string (Prim/Fwd/Rev/etc.)
        :param name: C++ name of variable
        :param size: Number of simplex elements
        :param value: Scalar value to initialize matrix with
        """
        self.overload = overload
        self.name = name
        self.stan_arg = "vector"
        self.size = size

        if value is None:
            self.value = 0.4
        else:
            self.value = value
    
    def is_reverse_mode(self):
        """Return true if the overload is reverse mode"""
        return self.overload.startswith("Rev")

    def is_eigen_compatible(self):
        """Return true (simplices are vectors)"""
        return True
    
    def is_varmat_compatible(self):
        """Return true (vectors are varmat compatible)"""
        return True
    
    def cpp(self):
        """Generate C++"""
        scalar = overload_scalar[self.overload]
        cpp_arg_template = get_cpp_type(self.stan_arg)
        arg_type = cpp_arg_template.replace("SCALAR", scalar)

        return "auto {name} = stan::test::make_simplex<{type}>({value}, {size});".format(name = self.name, type = arg_type, value = self.value, size = self.size)

class PositiveDefiniteMatrixVariable(CppStatement):
    """Represents a positive definite matrix variable"""
    def __init__(self, overload, name, size, value = None):
        """
        Initialize positive definite matrix
        
        :param overload: Type of overload as string (Prim/Fwd/Rev/etc.)
        :param name: C++ name of variable
        :param size: Number of rows/columns of positive definite matrix
        :param value: Scalar value to initialize matrix with
        """
        self.overload = overload
        self.name = name
        self.stan_arg = "matrix"
        self.size = size

        if value is None:
            self.value = 0.4
        else:
            self.value = value
    
    def is_reverse_mode(self):
        """Return true if the overload is reverse mode"""
        return self.overload.startswith("Rev")

    def is_eigen_compatible(self):
        """Return true (positive definite matrices are matrices)"""
        return True
    
    def is_varmat_compatible(self):
        """Return true (positive definite matrices are varmat compatible)"""
        return True
    
    def cpp(self):
        """Generate C++"""
        scalar = overload_scalar[self.overload]
        cpp_arg_template = get_cpp_type(self.stan_arg)
        arg_type = cpp_arg_template.replace("SCALAR", scalar)

        return "auto {name} = stan::test::make_pos_definite_matrix<{type}>({value}, {size});".format(name = self.name, type = arg_type, value = self.value, size = self.size)

class AlgebraSolverFunctorVariable(CppStatement):
    """Represents a functor variable that is compatible with Stan's algebra solver"""
    def __init__(self, name):
        """
        Initialize functor for algebra solver
        
        :param overload: Type of overload as string (Prim/Fwd/Rev/etc.)
        :param name: C++ name of variable
        """
        self.name = name
    
    def cpp(self):
        """Generate C++"""
        return "stan::test::simple_eq_functor {name};".format(name = self.name)

class OdeFunctorVariable(CppStatement):
    """Represents a functor variable that is compatible with Stan's variadic ode solvers"""
    def __init__(self, name):
        """
        Initialize functor for ode function
        
        :param overload: Type of overload as string (Prim/Fwd/Rev/etc.)
        :param name: C++ name of variable
        """
        self.name = name
    
    def cpp(self):
        """Generate C++"""
        return "stan::test::test_functor {name};".format(name = self.name)

class ReturnTypeTVariable(CppStatement):
    """Represents a scalar variable with type computed from the types of all the input arguments and stan::return_type_t"""
    def __init__(self, name, *args):
        """
        Initialize a scalar return type variable with type computed from the given args
        
        :param name: C++ name of variable
        :param args: Args from which to compute the return type
        """
        self.name = name
        arg_names = ["decltype({name})".format(name = arg.name) for arg in args]
        self.type_str = "stan::return_type_t<{names}>".format(names = ','.join(arg_names))
        self._is_reverse_mode = any(arg.is_reverse_mode() for arg in args)
    
    def is_reverse_mode(self):
        """Return true if any of the arguments are reverse mode"""
        return self._is_reverse_mode
    
    def cpp(self):
        """Generate C++"""
        return "{type} {name} = 0;".format(type = self.type_str, name = self.name)

class RngVariable(CppStatement):
    """Represents a random number generator"""
    def __init__(self, name):
        """
        Initialize an rng object
        
        :param name: C++ name of variable
        """
        self.name = name
    
    def cpp(self):
        """Generate C++"""
        return "std::minstd_rand {name};".format(name = self.name)

class OStreamVariable(CppStatement):
    """Represents an ostream pointer"""
    def __init__(self, name):
        """
        Initialize an ostream pointer
        
        :param name: C++ name of variable
        """
        self.name = name
    
    def cpp(self):
        """Generate C++"""
        return "std::ostream* {name} = &std::cout;".format(name = self.name)

class ArrayVariable(CppStatement):
    """Represents an array variable"""
    def __init__(self, overload, name, number_nested_arrays, inner_type, size, value):
        """
        Initialize an array
        
        :param overload: Type of overload as string (Prim/Fwd/Rev/etc.)
        :param name: C++ name of variable
        :param number_nested_arrays: Number of array dimensions
        :param inner_type: Stanc3 type string of inner type of array
        :param size: Size of arrays
        :param value: Either a scalar value (represented as another CppStatement) to initialize the array with or an iterable type of Python values (not CppStatements)
        """
        self.overload = overload
        self.name = name
        self.number_nested_arrays = number_nested_arrays
        self.inner_type = inner_type
        self.value = value
        self.size = size
    
    def is_reverse_mode(self):
        """
        If the underlying value is a CppStatement, return the reverse-mode-ness of that underlying variable

        Otherwise return true if the overload is reverse mode and inner type isn't an integer
        """
        if isinstance(self.value, CppStatement):
            return self.value.is_reverse_mode()
        else:
            return self.overload.startswith("Rev") and self.inner_type != "int"

    def is_varmat_compatible(self):
        """Return true if the underlying value is a CppStatement and is varmat compatible"""
        if isinstance(self.value, CppStatement):
            return self.value.is_varmat_compatible()
        else:
            return False
    
    def cpp(self):
        """Generate C++"""
        if isinstance(self.value, CppStatement):
            lhs_string = "decltype({name})".format(name = self.value.name)
            rhs_string = self.value.name
            for n in range(self.number_nested_arrays):
                lhs_string = "std::vector<" + lhs_string + ">"
                rhs_string = "{" + ",".join([rhs_string] * self.size) + "}"

            return "{lhs_string} {name} = {rhs_string};".format(lhs_string = lhs_string, name = self.name, rhs_string = rhs_string)
        else:
            # For when self.value is an initializer list
            if self.number_nested_arrays == 1:
                scalar = overload_scalar[self.overload]
                lhs_string = "std::vector<{type}>".format(type = arg_types[self.inner_type].replace("SCALAR", scalar))

                value_string = ",".join([repr(value) for value in self.value])
                return "{lhs_string} {name} = {{{value}}};".format(lhs_string = lhs_string, name = self.name, value = value_string)
            else:
                raise NotImplementedError("Array initializers are not implemented for number_nested_arrays = " + repr(self.number_nested_arrays))

class FunctionCall(CppStatement):
    """Represents a function call and optional assignment"""
    def __init__(self, function_name, name, *args):
        """
        Represents a function call and optional assignment of the result to a variable
        
        :param function_name: C++ name of function
        :param name: C++ name of variable, None if result is not saved
        :param args: Args to pass to function call
        """
        self.function_name = function_name
        self.name = name
        self.arg_str = ','.join(arg.name for arg in args)
        self._is_reverse_mode = any(arg.is_reverse_mode() for arg in args) and self.name is not None

    def is_reverse_mode(self):
        """
        Return true if any of the inputs are reverse mode and the output is saved

        This isn't always accurate cause the output could be an int or a boolean
        or something even if the input is reverse mode
        """
        return self._is_reverse_mode

    def is_eigen_compatible(self):
        raise Exception("is_eigen_compatible not implemented for FunctionCallAssign")

    def is_varmat(self):
        raise Exception("is_varmat not implemented for FunctionCallAssign")

    def cpp(self):
        """Generate C++"""
        if self.name:
            return "auto {name} = stan::math::eval({function_name}({arg_str}));".format(name = self.name, function_name = self.function_name, arg_str = self.arg_str)
        else:
            return "{function_name}({arg_str});".format(function_name = self.function_name, arg_str = self.arg_str)

class ExpressionVariable(CppStatement):
    """Represents an Eigen expression"""
    def __init__(self, name, arg, size = None):
        """
        Initialize Eigen expression variable
        
        :param name: C++ name of variable
        :param arg: Variable from which to form the expression
        :param size: Length of vector expression or rows/columns of matrix expression, if None use arg.size
        """
        self.name = name
        self.counter = IntVariable("{name}_counter".format(name = self.name), 0)
        if not size:
            self.size = arg.size
        else:
            self.size = size
        self._arg_overload = arg.overload
        self._is_reverse_mode = arg.is_reverse_mode()
        self._arg_stan_arg = arg.stan_arg
        self._arg_name = arg.name
    
    def is_reverse_mode(self):
        """Return true if the base argument was reverse mode"""
        return self._is_reverse_mode

    def is_eigen_compatible(self):
        """Return true (expressions are Eigen types)"""
        return True
    
    def is_expression(self):
        """Return true (this is an expression)"""
        return True

    def cpp(self):
        """Generate C++"""

        scalar = overload_scalar[self._arg_overload]
        counter_op_name = "{name}_counter_op".format(name = self.name)

        code = (
            self.counter.cpp() + os.linesep +
            "stan::test::counterOp<{scalar}> {counter_op_name}(&{counter});".format(scalar = scalar, counter_op_name = counter_op_name, counter = self.counter.name) + os.linesep
        )

        if self._arg_stan_arg == "matrix":
            return code + "auto {name} = {arg}.block(0,0,{size},{size}).unaryExpr({counter_op});".format(name = self.name, arg = self._arg_name, size = self.size, counter_op = counter_op_name)
        elif self._arg_stan_arg in ("vector", "row_vector"):
            return code + "auto {name} = {arg}.segment(0,{size}).unaryExpr({counter_op});".format(name = self.name, arg = self._arg_name, size = self.size, counter_op = counter_op_name)
        else:
            raise Exception("Can't make an expression out of a " + self.arg.stan_arg)
