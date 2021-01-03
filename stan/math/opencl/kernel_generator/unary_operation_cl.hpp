#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_UNARY_OPERATION_CL_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_UNARY_OPERATION_CL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <string>
#include <type_traits>
#include <set>
#include <utility>

namespace stan {
namespace math {

/**
 * Represents a unary operation in kernel generator expressions.
 * @tparam Derived derived type
 * @tparam T type of argument
 */
template <typename Derived, typename T, typename Scal>
class unary_operation_cl
    : public operation_cl<Derived, typename std::remove_reference_t<T>::Scalar,
                          T> {
 public:
  using Scalar = Scal;
  using base = operation_cl<Derived, Scalar, T>;
  using base::var_name_;

  /**
   * Constructor
   * @param a argument expression
   * @param op operation
   */
  unary_operation_cl(T&& a, const std::string& op)
      : base(std::forward<T>(a)), op_(op) {}

  /**
   * Generates kernel code for this expression.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether whether caller already handled matrix view
   * @param var_name_arg variable name of the nested expression
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& row_index_name,
                               const std::string& col_index_name,
                               const bool view_handled,
                               const std::string& var_name_arg) const {
    kernel_parts res{};
    res.body = type_str<Scalar>() + " " + var_name_ + " = " + op_ + var_name_arg
               + ";\n";
    return res;
  }

 protected:
  std::string op_;
};

/**
 * Represents a logical negation in kernel generator expressions.
 * @tparam Derived derived type
 * @tparam T type of argument
 */
template <typename T>
class logical_negation_
    : public unary_operation_cl<logical_negation_<T>, T, bool> {
  static_assert(
      std::is_integral<typename std::remove_reference_t<T>::Scalar>::value,
      "logical_negation: argument must be expression with integral "
      "or boolean return type!");
  using base = unary_operation_cl<logical_negation_<T>, T, bool>;
  using base::arguments_;

 public:
  /**
   * Constructor
   * @param a argument expression
   */
  explicit logical_negation_(T&& a) : base(std::forward<T>(a), "!") {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return logical_negation_<std::remove_reference_t<decltype(arg_copy)>>{
        std::move(arg_copy)};
  }

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const { return matrix_cl_view::Entire; }
};

/**
 * Logical negation of a kernel generator expression.
 *
 * @tparam T type of the argument
 * @param a argument expression
 * @return logical negation of given expression
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline logical_negation_<as_operation_cl_t<T>> operator!(T&& a) {
  return logical_negation_<as_operation_cl_t<T>>(
      as_operation_cl(std::forward<T>(a)));
}

/**
 * Represents an unary minus operation in kernel generator expressions.
 * @tparam Derived derived type
 * @tparam T type of argument
 */
template <typename T>
class unary_minus_
    : public unary_operation_cl<unary_minus_<T>, T,
                                typename std::remove_reference_t<T>::Scalar> {
  using base = unary_operation_cl<unary_minus_<T>, T,
                                  typename std::remove_reference_t<T>::Scalar>;
  using base::arguments_;

 public:
  /**
   * Constructor
   * @param a argument expression
   */
  explicit unary_minus_(T&& a) : base(std::forward<T>(a), "-") {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return unary_minus_<std::remove_reference_t<decltype(arg_copy)>>{
        std::move(arg_copy)};
  }

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const {
    return this->template get_arg<0>().view();
  }
};

/**
 * Unary minus of a kernel generator expression.
 *
 * @tparam T type of the argument
 * @param a argument expression
 * @return unary minus of given expression
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline unary_minus_<as_operation_cl_t<T>> operator-(T&& a) {
  return unary_minus_<as_operation_cl_t<T>>(
      as_operation_cl(std::forward<T>(a)));
}

}  // namespace math
}  // namespace stan

#endif
#endif
