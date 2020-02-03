#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_BINARY_OPERATOR_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_BINARY_OPERATOR_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/err.hpp>
#include <stan/math/opencl/multiply.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/scalar.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <stan/math/opencl/kernel_generator/common_return_scalar.hpp>
#include <string>
#include <type_traits>
#include <set>
#include <utility>

namespace stan {
namespace math {

/**
 * Represents a binary operation in kernel generator expressions.
 * @tparam Derived derived type
 * @tparam T_res scalar type of the result
 * @tparam T_a type of first argument
 * @tparam T_b type of second argument
 */
template <typename Derived, typename T_res, typename T_a, typename T_b>
class binary_operation : public operation_cl<Derived, T_res, T_a, T_b> {
 public:
  using Scalar = T_res;
  using base = operation_cl<Derived, Scalar, T_a, T_b>;
  using base::var_name;

 protected:
  std::string op_;
  using base::arguments_;

 public:
  /**
   * Constructor
   * @param a first argument of the binary operation
   * @param b sedond argument of the binary operation
   * @param op operator used in the binary operation
   */
  binary_operation(T_a&& a, T_b&& b, const std::string& op)  // NOLINT
      : base(std::forward<T_a>(a), std::forward<T_b>(b)), op_(op) {
    const std::string function = "binary_operator" + op;
    if (a.rows() != base::dynamic && b.rows() != base::dynamic) {
      check_size_match(function.c_str(), "Rows of ", "a", a.rows(), "rows of ",
                       "b", b.rows());
    }
    if (a.cols() != base::dynamic && b.cols() != base::dynamic) {
      check_size_match(function.c_str(), "Columns of ", "a", a.cols(),
                       "columns of ", "b", b.cols());
    }
  }

  /**
   * generates kernel code for this expression.
   * @param i row index variable name
   * @param j column index variable name
   * @param var_name_a variable name of the first nested expression
   * @param var_name_b variable name of the second nested expression
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& i, const std::string& j,
                               const std::string& var_name_a,
                               const std::string& var_name_b) const {
    kernel_parts res{};
    res.body = type_str<Scalar>() + " " + var_name + " = " + var_name_a + " "
               + op_ + " " + var_name_b + ";\n";
    return res;
  }

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const {
    return either(std::get<0>(arguments_).view(),
                  std::get<1>(arguments_).view());
  }
};

/**
  This macro is used to allow passing comma as part of a parameter in another
  macro.
  */
#define COMMA ,

/**
  Defines a new binary operation in kernel generator.
  @param class_name The name of the class this macro will define to represent
  this operation
  @param function_name The name of the function this macro will define that will
  be used to create this operation.
  @param scalar_type_expr The type of the scalar in the result of this
  operation. Can be a C++ expression that uses \c T_a and \c T_b as types of the
  scalars in the arguments to this operation.
  @param operation String containing operator that is used to implement this
  operation in kernel. Should be a valid infix operator in OpenCL C.
  */
#define ADD_BINARY_OPERATION(class_name, function_name, scalar_type_expr,     \
                             operation)                                       \
  template <typename T_a, typename T_b>                                       \
  class class_name : public binary_operation<class_name<T_a, T_b>,            \
                                             scalar_type_expr, T_a, T_b> {    \
   public:                                                                    \
    class_name(T_a&& a, T_b&& b) /* NOLINT */                                 \
        : binary_operation<class_name<T_a, T_b>, scalar_type_expr, T_a, T_b>( \
              std::forward<T_a>(a), std::forward<T_b>(b), operation) {}       \
  };                                                                          \
                                                                              \
  template <typename T_a, typename T_b,                                       \
            typename = require_all_valid_expressions_t<T_a, T_b>>             \
  inline class_name<as_operation_cl_t<T_a>, as_operation_cl_t<T_b>>           \
  function_name(T_a&& a, T_b&& b) { /* NOLINT */                              \
    return {as_operation_cl(std::forward<T_a>(a)),                            \
            as_operation_cl(std::forward<T_b>(b))};                           \
  }

/**
  Defines a new binary operation in kernel generator that needs to implement
  custom function that determines the view of the result.
  @param class_name The name of the class this macro will define to represent
  this operation
  @param function_name The name of the function this macro will define that will
  be used to create this operation.
  @param scalar_type_expr The type of the scalar in the result of this
  operation. Can be a C++ expression that uses \c T_a and \c T_b as types of the
  scalars in the arguments to this operation.
  @param operation String containing operator that is used to implement this
  operation in kernel. Should be a valid infix operator in OpenCL C.
  @param ... Code that implements body of the \c .view() member function of the
  class that represents this expression. Should return an object of type
  matrix_cl_view. Can use \c base::arguments_ to access arguments to this
  expression. This is a variadic argument to allow commas in code with no
  special handling.
  */
#define ADD_BINARY_OPERATION_WITH_CUSTOM_VIEW(                                \
    class_name, function_name, scalar_type_expr, operation, ...)              \
  template <typename T_a, typename T_b>                                       \
  class class_name : public binary_operation<class_name<T_a, T_b>,            \
                                             scalar_type_expr, T_a, T_b> {    \
   public:                                                                    \
    class_name(T_a&& a, T_b&& b) /* NOLINT */                                 \
        : binary_operation<class_name<T_a, T_b>, scalar_type_expr, T_a, T_b>( \
              std::forward<T_a>(a), std::forward<T_b>(b), operation) {}       \
    inline matrix_cl_view view() const { __VA_ARGS__; }                       \
  };                                                                          \
                                                                              \
  template <typename T_a, typename T_b,                                       \
            typename = require_all_valid_expressions_t<T_a, T_b>>             \
  inline class_name<as_operation_cl_t<T_a>, as_operation_cl_t<T_b>>           \
  function_name(T_a&& a, T_b&& b) { /* NOLINT */                              \
    return {as_operation_cl(std::forward<T_a>(a)),                            \
            as_operation_cl(std::forward<T_b>(b))};                           \
  }

ADD_BINARY_OPERATION(addition_, operator+, common_scalar_t<T_a COMMA T_b>, "+");
ADD_BINARY_OPERATION(subtraction_, operator-, common_scalar_t<T_a COMMA T_b>,
                     "-");
ADD_BINARY_OPERATION_WITH_CUSTOM_VIEW(
    elewise_multiplication_, elewise_multiplication,
    common_scalar_t<T_a COMMA T_b>, "*",
    using base = binary_operation<elewise_multiplication_<T_a, T_b>,
                                  common_scalar_t<T_a, T_b>, T_a, T_b>;
    return both(std::get<0>(base::arguments_).view(),
                std::get<1>(base::arguments_).view()););
ADD_BINARY_OPERATION_WITH_CUSTOM_VIEW(
    elewise_division_, elewise_division, common_scalar_t<T_a COMMA T_b>, "/",
    using base = binary_operation<elewise_division_<T_a, T_b>,
                                  common_scalar_t<T_a, T_b>, T_a, T_b>;
    return either(std::get<0>(base::arguments_).view(),
                  invert(std::get<1>(base::arguments_).view())););
ADD_BINARY_OPERATION(less_than_, operator<, bool, "<");
ADD_BINARY_OPERATION_WITH_CUSTOM_VIEW(less_than_or_equal_, operator<=, bool,
                                      "<=", return matrix_cl_view::Entire);
ADD_BINARY_OPERATION(greater_than_, operator>, bool, ">");
ADD_BINARY_OPERATION_WITH_CUSTOM_VIEW(greater_than_or_equal_, operator>=, bool,
                                      ">=", return matrix_cl_view::Entire);
ADD_BINARY_OPERATION_WITH_CUSTOM_VIEW(equals_, operator==, bool,
                                      "==", return matrix_cl_view::Entire);
ADD_BINARY_OPERATION(not_equals_, operator!=, bool, "!=");

ADD_BINARY_OPERATION(logical_or_, operator||, bool, "||");
ADD_BINARY_OPERATION(logical_and_, operator&&, bool, "&&");

/**
 * Multiplication of a scalar and a kernel generator expression.
 * @tparam T_a type of scalar
 * @tparam T_b type of expression
 * @param a scalar
 * @param b expression
 * @return Multiplication of given arguments
 */
template <typename T_a, typename T_b, typename = require_arithmetic_t<T_a>,
          typename = require_all_valid_expressions_t<T_b>>
inline elewise_multiplication_<scalar_<T_a>, as_operation_cl_t<T_b>> operator*(
    T_a&& a, T_b&& b) {  // NOLINT
  return {as_operation_cl(std::forward<T_a>(a)),
          as_operation_cl(std::forward<T_b>(b))};
}

/**
 * Multiplication of a kernel generator expression and a scalar.
 * @tparam T_a type of expression
 * @tparam T_b type of scalar
 * @param a expression
 * @param b scalar
 * @return Multiplication of given arguments
 */
template <typename T_a, typename T_b,
          typename = require_all_valid_expressions_t<T_a>,
          typename = require_arithmetic_t<T_b>>
inline elewise_multiplication_<as_operation_cl_t<T_a>, scalar_<T_b>> operator*(
    T_a&& a, const T_b b) {  // NOLINT
  return {as_operation_cl(std::forward<T_a>(a)), as_operation_cl(b)};
}

/**
 * Matrix multiplication of two kernel generator expressions. Evaluates both
 * expressions before calculating the matrix product.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return Matrix product of given arguments
 */
template <typename T_a, typename T_b,
          typename = require_all_valid_expressions_and_none_scalar_t<T_a, T_b>>
inline matrix_cl<double> operator*(const T_a& a, const T_b& b) {
  // no need for perfect forwarding as operations are evaluated
  return stan::math::opencl::multiply(as_operation_cl(a).eval(),
                                      as_operation_cl(b).eval());
}

#undef COMMA
#undef ADD_BINARY_OPERATION
#undef ADD_BINARY_OPERATION_WITH_CUSTOM_VIEW

}  // namespace math
}  // namespace stan
#endif
#endif
