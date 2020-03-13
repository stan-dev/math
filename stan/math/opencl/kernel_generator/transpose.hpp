#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_TRANSPOSE_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_TRANSPOSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <algorithm>
#include <string>
#include <set>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {

/**
 * Represents a transpose in kernel generator expressions.
 *
 * @tparam Derived derived type
 * @tparam Arg type of the argument
 */
template <typename Arg>
class transpose_
    : public operation_cl<transpose_<Arg>,
                          typename std::remove_reference_t<Arg>::Scalar, Arg> {
 public:
  using Scalar = typename std::remove_reference_t<Arg>::Scalar;
  using base = operation_cl<transpose_<Arg>, Scalar, Arg>;
  using base::var_name;

 protected:
  using base::arguments_;

 public:
  /**
   * Constructor
   * @param a expression to transpose
   */
  explicit transpose_(Arg&& a) : base(std::forward<Arg>(a)) {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline transpose_<std::remove_reference_t<Arg>> deep_copy() {
    return transpose_<std::remove_reference_t<Arg>>{
        std::get<0>(arguments_).deep_copy()};
  }

  /**
   * generates kernel code for this and nested expressions.
   * @param[in,out] generated set of already generated operations
   * @param ng name generator for this kernel
   * @param i row index variable name
   * @param j column index variable name
   * @return part of kernel with code for this and nested expressions
   */
  inline kernel_parts generate(const std::string& i, const std::string& j,
                               const std::string& var_name_arg) const {
    var_name = var_name_arg;
    return {};
  }

  /**
   * Swaps indices \c i and \c j for the argument expression.
   * @param[in, out] i row index
   * @param[in, out] j column index
   */
  inline void modify_argument_indices(std::string& i, std::string& j) const {
    std::swap(i, j);
  }

  /**
   * Number of rows of a matrix that would be the result of evaluating this
   * expression.
   * @return number of rows
   */
  inline int rows() const { return std::get<0>(arguments_).cols(); }

  /**
   * Number of columns of a matrix that would be the result of evaluating this
   * expression.
   * @return number of columns
   */
  inline int cols() const { return std::get<0>(arguments_).rows(); }

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const {
    return transpose(std::get<0>(arguments_).view());
  }

  /**
   * Determine index of bottom diagonal written.
   * @return index of bottom diagonal
   */
  inline int bottom_diagonal() const {
    return std::get<0>(arguments_).top_diagonal();
  }

  /**
   * Determine index of top diagonal written.
   * @return index of top diagonal
   */
  inline int top_diagonal() const {
    return std::get<0>(arguments_).bottom_diagonal();
  }
};

template <typename Arg,
          typename = require_all_valid_expressions_and_none_scalar_t<Arg>>
inline auto transpose(Arg&& a) {
  auto&& a_operation = as_operation_cl(std::forward<Arg>(a)).deep_copy();
  return transpose_<std::remove_reference_t<decltype(a_operation)>>{
      std::move(a_operation)};
}

}  // namespace math
}  // namespace stan

#endif
#endif
