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
 * @tparam T_a type of first argument
 * @tparam T_b type of second argument
 */
template <typename Derived, typename T_a, typename T_b>
class binary_operation
    : public operation_cl<Derived, common_scalar_t<T_a, T_b>, T_a, T_b> {
 public:
  using Scalar = common_scalar_t<T_a, T_b>;
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
 * Represents addition in kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 */
template <typename T_a, typename T_b>
class addition_ : public binary_operation<addition_<T_a, T_b>, T_a, T_b> {
 public:
  /**
   * Constructor.
   * @param a first expression
   * @param b second expression
   */
  addition_(T_a&& a, T_b&& b)  // NOLINT
      : binary_operation<addition_<T_a, T_b>, T_a, T_b>(
            std::forward<T_a>(a), std::forward<T_b>(b), "+") {}
};

/**
 * Addition of two kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first argument
 * @param b second argument
 * @return Addition of given expressions
 */
template <typename T_a, typename T_b,
          typename = require_all_valid_expressions_t<T_a, T_b>>
inline addition_<as_operation_cl_t<T_a>, as_operation_cl_t<T_b>> operator+(
    T_a&& a, T_b&& b) {  // NOLINT
  return {as_operation_cl(std::forward<T_a>(a)),
          as_operation_cl(std::forward<T_b>(b))};
}

/**
 * Represents subtraction in kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 */
template <typename T_a, typename T_b>
class subtraction_ : public binary_operation<subtraction_<T_a, T_b>, T_a, T_b> {
 public:
  /**
   * Constructor.
   * @param a first expression
   * @param b second expression
   */
  subtraction_(T_a&& a, T_b&& b)  // NOLINT
      : binary_operation<subtraction_<T_a, T_b>, T_a, T_b>(
            std::forward<T_a>(a), std::forward<T_b>(b), "-") {}
};

/**
 * Subtraction of two kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return Subtraction of given expressions
 */
template <typename T_a, typename T_b,
          typename = require_all_valid_expressions_t<T_a, T_b>>
inline subtraction_<as_operation_cl_t<T_a>, as_operation_cl_t<T_b>> operator-(
    T_a&& a, T_b&& b) {  // NOLINT
  return {as_operation_cl(std::forward<T_a>(a)),
          as_operation_cl(std::forward<T_b>(b))};
}

/**
 * Represents element-wise multiplication in kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 */
template <typename T_a, typename T_b,
          typename = require_all_valid_expressions_t<T_a, T_b>>
class elewise_multiplication_
    : public binary_operation<elewise_multiplication_<T_a, T_b>, T_a, T_b> {
 public:
  /**
   * Constructor.
   * @param a first expression
   * @param b second expression
   */
  elewise_multiplication_(T_a&& a, T_b&& b)  // NOLINT
      : binary_operation<elewise_multiplication_<T_a, T_b>, T_a, T_b>(
            std::forward<T_a>(a), std::forward<T_b>(b), "*") {}

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const {
    using base = binary_operation<elewise_multiplication_<T_a, T_b>, T_a, T_b>;
    return both(std::get<0>(base::arguments_).view(),
                std::get<1>(base::arguments_).view());
  }
};

/**
 * Element-wise multiplication of two kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return Element-wise multiplication of given expressions
 */
template <typename T_a, typename T_b>
inline elewise_multiplication_<as_operation_cl_t<T_a>, as_operation_cl_t<T_b>>
elewise_multiplication(T_a&& a, T_b&& b) {  // NOLINT
  return {as_operation_cl(std::forward<T_a>(a)),
          as_operation_cl(std::forward<T_b>(b))};
}

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

/**
 * Represents element-wise division in kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 */
template <typename T_a, typename T_b>
class elewise_division_
    : public binary_operation<elewise_division_<T_a, T_b>, T_a, T_b> {
 public:
  /**
   * Constructor.
   * @param a first expression
   * @param b second expression
   */
  elewise_division_(T_a&& a, T_b&& b)  // NOLINT
      : binary_operation<elewise_division_<T_a, T_b>, T_a, T_b>(
            std::forward<T_a>(a), std::forward<T_b>(b), "/") {}

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const {
    using base = binary_operation<elewise_division_<T_a, T_b>, T_a, T_b>;
    return either(std::get<0>(base::arguments_).view(),
                  invert(std::get<1>(base::arguments_).view()));
  }
};

/**
 * Element-wise division of two kernel generator expressions.
 * @tparam T_a type of first expression
 * @tparam T_b type of second expression
 * @param a first expression
 * @param b second expression
 * @return Element-wise division of given expressions
 */
template <typename T_a, typename T_b,
          typename = require_all_valid_expressions_t<T_a, T_b>>
inline elewise_division_<as_operation_cl_t<T_a>, as_operation_cl_t<T_b>>
elewise_division(T_a&& a, T_b&& b) {  // NOLINT
  return {as_operation_cl(std::forward<T_a>(a)),
          as_operation_cl(std::forward<T_b>(b))};
}

}  // namespace math
}  // namespace stan
#endif
#endif
