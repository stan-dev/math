#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_TRANSPOSE_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_TRANSPOSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl_lhs.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <algorithm>
#include <string>
#include <set>
#include <tuple>
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
    : public operation_cl_lhs<
          transpose_<Arg>, typename std::remove_reference_t<Arg>::Scalar, Arg> {
 public:
  using Scalar = typename std::remove_reference_t<Arg>::Scalar;
  using base = operation_cl_lhs<transpose_<Arg>, Scalar, Arg>;
  using base::var_name_;
  using view_transitivity = std::tuple<std::true_type>;
  using base::operator=;

  /**
   * Constructor
   * @param a expression to transpose
   */
  explicit transpose_(Arg&& a) : base(std::forward<Arg>(a)) {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return transpose_<std::remove_reference_t<decltype(arg_copy)>>{
        std::move(arg_copy)};
  }

  /**
   * Swaps indices \c row_index_name and \c col_index_name for the argument
   * expression.
   * @param[in, out] row_index_name row index
   * @param[in, out] col_index_name column index
   */
  inline void modify_argument_indices(std::string& row_index_name,
                                      std::string& col_index_name) const {
    std::swap(row_index_name, col_index_name);
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

  /**
   * Sets the view of the underlying matrix depending on which of its parts are
   * written to.
   * @param bottom_diagonal Index of the top sub- or super- diagonal written
   * with nonzero elements.
   * @param top_diagonal Index of the top sub- or super- diagonal written with
   * nonzero elements.
   * @param bottom_zero_diagonal Index of the top sub- or super- diagonal
   * written with zeros if it ie more extreme than \c bottom_diagonal. Otherwise
   * it should be set to equal value as \c bottom_diagonal.
   * @param top_zero_diagonal Index of the top sub- or super- diagonal written
   * with zeros if it ie more extreme than \c top_diagonal. Otherwise it should
   * be set to equal value as \c top_diagonal.
   */
  void set_view(int bottom_diagonal, int top_diagonal, int bottom_zero_diagonal,
                int top_zero_diagonal) const {
    this->template get_arg<0>().set_view(
        top_diagonal, bottom_diagonal, top_zero_diagonal, bottom_zero_diagonal);
  }

  /**
   * Sets the dimensions of the underlying expressions if possible. If not
   * checks whether they have correct dimensions.
   * @param rows desired number of rows
   * @param cols desired number of columns
   * @throws std::invalid_argument desired dimensions do not match with
   * dimensions of underlying expression that can not be resized.
   */
  void check_assign_dimensions(int rows, int cols) const {
    this->template get_arg<0>().check_assign_dimensions(cols, rows);
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
          typename = require_all_kernel_expressions_and_none_scalar_t<Arg>>
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
