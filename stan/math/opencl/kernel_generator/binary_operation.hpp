#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_BINARY_OPERATOR_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_BINARY_OPERATOR_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/scalar.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/common_return_scalar.hpp>
#include <algorithm>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */

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
  using base::var_name_;

 protected:
  std::string op_;
  using base::arguments_;

 public:
  /**
   * Constructor
   * @param a first argument of the binary operation
   * @param b second argument of the binary operation
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
   * Generates kernel code for this expression.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether whether caller already handled matrix view
   * @param var_name_a variable name of the first nested expression
   * @param var_name_b variable name of the second nested expression
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& row_index_name,
                               const std::string& col_index_name,
                               const bool view_handled,
                               const std::string& var_name_a,
                               const std::string& var_name_b) const {
    kernel_parts res{};
    res.body = type_str<Scalar>() + " " + var_name_ + " = " + var_name_a + " "
               + op_ + " " + var_name_b + ";\n";
    return res;
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
    using base                                                                \
        = binary_operation<class_name<T_a, T_b>, scalar_type_expr, T_a, T_b>; \
    using base::arguments_;                                                   \
                                                                              \
   public:                                                                    \
    using base::rows;                                                         \
    using base::cols;                                                         \
    class_name(T_a&& a, T_b&& b) /* NOLINT */                                 \
        : base(std::forward<T_a>(a), std::forward<T_b>(b), operation) {}      \
    inline auto deep_copy() const {                                           \
      auto&& a_copy = this->template get_arg<0>().deep_copy();                \
      auto&& b_copy = this->template get_arg<1>().deep_copy();                \
      return class_name<std::remove_reference_t<decltype(a_copy)>,            \
                        std::remove_reference_t<decltype(b_copy)>>(           \
          std::move(a_copy), std::move(b_copy));                              \
    }                                                                         \
  };                                                                          \
                                                                              \
  template <typename T_a, typename T_b,                                       \
            require_all_kernel_expressions_t<T_a, T_b>* = nullptr,            \
            require_any_not_arithmetic_t<T_a, T_b>* = nullptr>                \
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
#define ADD_BINARY_OPERATION_WITH_CUSTOM_CODE(                                \
    class_name, function_name, scalar_type_expr, operation, ...)              \
  template <typename T_a, typename T_b>                                       \
  class class_name : public binary_operation<class_name<T_a, T_b>,            \
                                             scalar_type_expr, T_a, T_b> {    \
    using base                                                                \
        = binary_operation<class_name<T_a, T_b>, scalar_type_expr, T_a, T_b>; \
    using base::arguments_;                                                   \
                                                                              \
   public:                                                                    \
    using base::rows;                                                         \
    using base::cols;                                                         \
    class_name(T_a&& a, T_b&& b) /* NOLINT */                                 \
        : base(std::forward<T_a>(a), std::forward<T_b>(b), operation) {}      \
    inline auto deep_copy() const {                                           \
      auto&& a_copy = this->template get_arg<0>().deep_copy();                \
      auto&& b_copy = this->template get_arg<1>().deep_copy();                \
      return class_name<std::remove_reference_t<decltype(a_copy)>,            \
                        std::remove_reference_t<decltype(b_copy)>>(           \
          std::move(a_copy), std::move(b_copy));                              \
    }                                                                         \
    __VA_ARGS__                                                               \
  };                                                                          \
                                                                              \
  template <typename T_a, typename T_b,                                       \
            require_all_kernel_expressions_t<T_a, T_b>* = nullptr,            \
            require_any_not_arithmetic_t<T_a, T_b>* = nullptr>                \
  inline class_name<as_operation_cl_t<T_a>, as_operation_cl_t<T_b>>           \
  function_name(T_a&& a, T_b&& b) { /* NOLINT */                              \
    return {as_operation_cl(std::forward<T_a>(a)),                            \
            as_operation_cl(std::forward<T_b>(b))};                           \
  }

ADD_BINARY_OPERATION(addition_operator_, operator+,
                     common_scalar_t<T_a COMMA T_b>, "+");
ADD_BINARY_OPERATION(addition_, add, common_scalar_t<T_a COMMA T_b>, "+");
ADD_BINARY_OPERATION(subtraction_operator_, operator-,
                     common_scalar_t<T_a COMMA T_b>, "-");
ADD_BINARY_OPERATION(subtraction_, subtract, common_scalar_t<T_a COMMA T_b>,
                     "-");
ADD_BINARY_OPERATION_WITH_CUSTOM_CODE(
    elt_multiply_, elt_multiply, common_scalar_t<T_a COMMA T_b>, "*",
    using view_transitivity = std::tuple<std::true_type, std::true_type>;
    inline std::pair<int, int> extreme_diagonals() const {
      std::pair<int, int> diags0
          = this->template get_arg<0>().extreme_diagonals();
      std::pair<int, int> diags1
          = this->template get_arg<1>().extreme_diagonals();
      return {std::max(diags0.first, diags1.first),
              std::min(diags0.second, diags1.second)};
    });
ADD_BINARY_OPERATION_WITH_CUSTOM_CODE(
    elt_divide_, elt_divide, common_scalar_t<T_a COMMA T_b>, "/",
    inline std::pair<int, int> extreme_diagonals() const {
      return {-rows() + 1, cols() - 1};
    });
ADD_BINARY_OPERATION_WITH_CUSTOM_CODE(
    elt_modulo_, operator%, common_scalar_t<T_a COMMA T_b>, "%",
    static_assert(
        std::is_integral<scalar_type_t<T_a>>::value&&
            std::is_integral<scalar_type_t<T_b>>::value,
        "both operands to operator% must have integral scalar types!");
    inline std::pair<int, int> extreme_diagonals() const {
      return {-rows() + 1, cols() - 1};
    });

ADD_BINARY_OPERATION(less_than_, operator<, bool, "<");
ADD_BINARY_OPERATION_WITH_CUSTOM_CODE(
    less_than_or_equal_, operator<=, bool,
    "<=", inline std::pair<int, int> extreme_diagonals() const {
      return {-rows() + 1, cols() - 1};
    });
ADD_BINARY_OPERATION(greater_than_, operator>, bool, ">");
ADD_BINARY_OPERATION_WITH_CUSTOM_CODE(
    greater_than_or_equal_, operator>=, bool,
    ">=", inline std::pair<int, int> extreme_diagonals() const {
      return {-rows() + 1, cols() - 1};
    });
ADD_BINARY_OPERATION_WITH_CUSTOM_CODE(
    equals_, operator==, bool,
    "==", inline std::pair<int, int> extreme_diagonals() const {
      return {-rows() + 1, cols() - 1};
    });
ADD_BINARY_OPERATION(not_equals_, operator!=, bool, "!=");

ADD_BINARY_OPERATION(logical_or_, operator||, bool, "||");
ADD_BINARY_OPERATION_WITH_CUSTOM_CODE(
    logical_and_, operator&&, bool, "&&",
    using view_transitivity = std::tuple<std::true_type, std::true_type>;
    inline std::pair<int, int> extreme_diagonals() const {
      std::pair<int, int> diags0
          = this->template get_arg<0>().extreme_diagonals();
      std::pair<int, int> diags1
          = this->template get_arg<1>().extreme_diagonals();
      return {std::max(diags0.first, diags1.first),
              std::min(diags0.second, diags1.second)};
    });

/**
 * Multiplication of a scalar and a kernel generator expression.
 * @tparam T_a type of scalar
 * @tparam T_b type of expression
 * @param a scalar
 * @param b expression
 * @return Multiplication of given arguments
 */
template <typename T_a, typename T_b, typename = require_arithmetic_t<T_a>,
          typename = require_all_kernel_expressions_t<T_b>>
inline elt_multiply_<as_operation_cl_t<T_a>, as_operation_cl_t<T_b>> operator*(
    T_a a, T_b&& b) {  // NOLINT
  return {as_operation_cl(a), as_operation_cl(std::forward<T_b>(b))};
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
          typename = require_all_kernel_expressions_t<T_a>,
          typename = require_arithmetic_t<T_b>>
inline elt_multiply_<as_operation_cl_t<T_a>, as_operation_cl_t<T_b>> operator*(
    T_a&& a, const T_b b) {  // NOLINT
  return {as_operation_cl(std::forward<T_a>(a)), as_operation_cl(b)};
}

#undef COMMA
#undef ADD_BINARY_OPERATION
#undef ADD_BINARY_OPERATION_WITH_CUSTOM_CODE
/** @}*/
}  // namespace math
}  // namespace stan
#endif
#endif
