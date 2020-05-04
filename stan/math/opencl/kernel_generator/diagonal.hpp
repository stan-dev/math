#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_DIAGONAL_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_DIAGONAL_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl_lhs.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <set>
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
  using base::var_name;
  using base::operator=;

  /**
   * Constructor
   * @param a expression
   */
  diagonal_(T&& a) : base(std::forward<T>(a)) {}

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
   * Generates kernel code for this and nested expressions.
   * @param[in,out] generated set of (pointer to) already generated operations
   * @param name_gen name generator for this kernel
   * @param i row index variable name
   * @param j column index variable name
   * @param view_handled whether caller already handled matrix view
   * @return part of kernel with code for this and nested expressions
   */
  inline kernel_parts get_kernel_parts(
      std::set<const operation_cl_base*>& generated, name_generator& name_gen,
      const std::string& i, const std::string& j, bool view_handled) const {
    kernel_parts res{};
    if (generated.count(this) == 0) {
      generated.insert(this);
      res = this->template get_arg<0>().get_kernel_parts(generated, name_gen, i,
                                                         i, true);
      var_name = this->template get_arg<0>().var_name;
    }
    return res;
  }

  /**
   * Sets j to value of i. This is only used when diagonal is assigned to.
   * @param[in, out] i row index
   * @param[in, out] j column index
   */
  inline void modify_argument_indices(std::string& i, std::string& j) const {
    j = i;
  }

  /**
   * Generates kernel code for this and nested expressions if this expression
   * appears on the left hand side of an assignment.
   * @param i row index variable name
   * @param j column index variable name
   * @param var_name_arg name of the variable in kernel that holds argument to
   * this expression
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate_lhs(const std::string& i, const std::string& j,
                                   const std::string& var_name_arg) const {
    return {};
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
   * Sets view of the underlying matrix depending on which part is written.
   * Setting view to diagonal never changes view of the underlying matrix.
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
          typename = require_all_valid_expressions_and_none_scalar_t<T>>
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
