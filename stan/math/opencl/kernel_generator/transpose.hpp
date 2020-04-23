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

/** \addtogroup opencl_kernel_generator
 *  @{
 */

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
  using view_transitivity = std::tuple<std::true_type>;

  /**
   * Constructor
   * @param a expression to transpose
   */
  explicit transpose_(Arg&& a) : base(std::forward<Arg>(a)) {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline transpose_<std::remove_reference_t<Arg>> deep_copy() const {
    return transpose_<std::remove_reference_t<Arg>>{
        this->template get_arg<0>().deep_copy()};
  }

  /**
   * generates kernel code for this and nested expressions.
   * @param i row index variable name.
   * @param j column index variable name.
   * @param view_handled whether whether caller already handled matrix view
   * @param var_name_arg The name of this variable.
   * @return part of kernel with code for this and nested expressions
   */
  inline kernel_parts generate(const std::string& i, const std::string& j,
                               const bool view_handeled,
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
  inline int rows() const { return this->template get_arg<0>().cols(); }

  /**
   * Number of columns of a matrix that would be the result of evaluating this
   * expression.
   * @return number of columns
   */
  inline int cols() const { return this->template get_arg<0>().rows(); }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    std::pair<int, int> arg_diags
        = this->template get_arg<0>().extreme_diagonals();
    return {-arg_diags.second, -arg_diags.first};
  }
};

/**
 * Transposes a kernel generator expression.
 *
 * Transposition modifies how its argument is indexed. If a matrix is both an
 * argument and result of such an operation (such as in <code> a = transpose(a);
 * </code>), the result will be wrong due to aliasing. In such case the
 * expression should be evaluating in a temporary by doing <code> a =
 * transpose(a).eval();</code>.
 * @tparam Arg type of the argument expression.
 * @param a argument to transposition
 */
template <typename Arg,
          typename = require_all_valid_expressions_and_none_scalar_t<Arg>>
inline auto transpose(Arg&& a) {
  auto&& a_operation = as_operation_cl(std::forward<Arg>(a)).deep_copy();
  return transpose_<std::remove_reference_t<decltype(a_operation)>>{
      std::move(a_operation)};
}
/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
