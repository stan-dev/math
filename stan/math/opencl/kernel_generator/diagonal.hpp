#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_DIAGONAL_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_DIAGONAL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl_lhs.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <algorithm>
#include <map>
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
 * Represents diagonal of a matrix (as column vector) in kernel generator
 * expressions.
 * @tparam Derived derived type
 * @tparam T type of argument
 */
template <typename T>
class diagonal_
    : public operation_cl_lhs<diagonal_<T>,
                              typename std::remove_reference_t<T>::Scalar, T> {
 public:
  using Scalar = typename std::remove_reference_t<T>::Scalar;
  using base = operation_cl_lhs<diagonal_<T>, Scalar, T>;
  using base::var_name_;
  using base::operator=;

  /**
   * Constructor
   * @param a expression
   */
  explicit diagonal_(T&& a) : base(std::forward<T>(a)) {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return diagonal_<std::remove_reference_t<decltype(arg_copy)>>{
        std::move(arg_copy)};
  }

  /**
   * Sets col_index_name to value of row_index_name. This is only used when
   * diagonal is assigned to.
   * @param[in, out] row_index_name row index
   * @param[in, out] col_index_name column index
   */
  inline void modify_argument_indices(std::string& row_index_name,
                                      std::string& col_index_name) const {
    col_index_name = row_index_name;
  }

  /**
   * Number of rows of a matrix that would be the result of evaluating this
   * expression.
   * @return number of rows
   */
  inline int rows() const {
    return std::min(this->template get_arg<0>().rows(),
                    this->template get_arg<0>().cols());
  }

  /**
   * Number of columns of a matrix that would be the result of evaluating this
   * expression.
   * @return 1
   */
  inline int cols() const { return 1; }

  /**
   * Sets the view of the underlying matrix depending on which of its parts are
   * written to. Setting view to diagonal never changes view of the underlying
   * matrix.
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
  inline void set_view(int bottom_diagonal, int top_diagonal,
                       int bottom_zero_diagonal, int top_zero_diagonal) const {}

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    return {1 - rows(), 1};
  }

  /**
   * Checks if desired dimensions match dimensions of the block.
   * @param rows desired number of rows
   * @param cols desired number of columns
   * @throws std::invalid_argument desired dimensions do not match dimensions
   * of the block.
   */
  inline void check_assign_dimensions(int rows, int cols) const {
    check_size_match("diagonal_.check_assign_dimensions", "Rows of ",
                     "diagonal", this->rows(), "rows of ", "expression", rows);
    check_size_match("diagonal_.check_assign_dimensions", "Columns of ",
                     "diagonal", 1, "columns of ", "expression", cols);
  }
};

/**
 * Diagonal of a kernel generator expression. Diagonal is represented as a
 * column vector.
 *
 * @tparam T type of argument
 * @param a input argument
 * @return Diagonal of given expression
 */
template <typename T,
          typename = require_all_kernel_expressions_and_none_scalar_t<T>>
inline auto diagonal(T&& a) {
  auto&& a_operation = as_operation_cl(std::forward<T>(a)).deep_copy();
  return diagonal_<std::remove_reference_t<decltype(a_operation)>>(
      std::move(a_operation));
}
/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
