#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_BLOCK_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_BLOCK_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl_lhs.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <set>
#include <string>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {

/**
 * Represents submatrix block in kernel generator expressions.
 * @tparam Derived derived type
 * @tparam T type of argument
 */
template <typename T>
class block_
    : public operation_cl_lhs<block_<T>,
                              typename std::remove_reference_t<T>::Scalar, T> {
 public:
  using Scalar = typename std::remove_reference_t<T>::Scalar;
  using base = operation_cl_lhs<block_<T>, Scalar, T>;
  using base::var_name;

 protected:
  int start_row_, start_col_, rows_, cols_;

 public:
  /**
   * Constructor
   * @param a expression
   * @param start_row first row of block
   * @param start_col first column of a block
   * @param rows number of rows in block
   * @param cols number of columns in block
   */
  block_(T&& a, int start_row, int start_col, int rows, int cols)
      : base(std::forward<T>(a)),
        start_row_(start_row),
        start_col_(start_col),
        rows_(rows),
        cols_(cols) {
    if ((a.rows() != base::dynamic && (start_row + rows) > a.rows())
        || (a.cols() != base::dynamic && (start_col + cols) > a.cols())) {
      throw_domain_error("block", "block of \"a\"", " is out of bounds", "");
    }
  }

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return block_<std::remove_reference_t<decltype(arg_copy)>>{
        std::move(arg_copy), start_row_, start_col_, rows_, cols_};
  }

  /**
   * Generates kernel code for this expression.
   * @param i row index variable name
   * @param j column index variable name
   * @param var_name_arg name of the variable in kernel that holds argument to
   * this expression
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& i, const std::string& j,
                               const std::string& var_name_arg) const {
    kernel_parts res;
    res.body
        = type_str<Scalar>() + " " + var_name + " = " + var_name_arg + ";\n";
    res.args = "int " + var_name + "_i, int " + var_name + "_j, ";
    return res;
  }

  /**
   * Sets offset of block to indices of the argument expression
   * @param[in, out] i row index
   * @param[in, out] j column index
   */
  inline void modify_argument_indices(std::string& i, std::string& j) const {
    i = "(" + i + " + " + var_name + "_i)";
    j = "(" + j + " + " + var_name + "_j)";
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
    kernel_parts res;
    res.args = "int " + var_name + "_i, int " + var_name + "_j, ";
    return res;
  }

  /**
   * Sets kernel arguments for this and nested expressions.
   * @param[in,out] generated set of expressions that already set their kernel
   * arguments
   * @param kernel kernel to set arguments on
   * @param[in,out] arg_num consecutive number of the first argument to set.
   * This is incremented for each argument set by this function.
   */
  inline void set_args(std::set<const operation_cl_base*>& generated,
                       cl::Kernel& kernel, int& arg_num) const {
    if (generated.count(this) == 0) {
      generated.insert(this);
      this->template get_arg<0>().set_args(generated, kernel, arg_num);
      kernel.setArg(arg_num++, start_row_);
      kernel.setArg(arg_num++, start_col_);
    }
  }

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const {
    matrix_cl_view view;
    if (bottom_diagonal() < 0) {
      view = matrix_cl_view::Lower;
    } else {
      view = matrix_cl_view::Diagonal;
    }
    if (top_diagonal() > 0) {
      view = either(view, matrix_cl_view::Upper);
    }
    return view;
  }

  /**
   * Number of rows of a matrix that would be the result of evaluating this
   * expression.
   * @return number of rows
   */
  inline int rows() const { return rows_; }

  /**
   * Number of columns of a matrix that would be the result of evaluating this
   * expression.
   * @return number of columns
   */
  inline int cols() const { return cols_; }

  /**
   * Sets view of the underlying matrix depending on which part is written.
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
                       int bottom_zero_diagonal, int top_zero_diagonal) const {
    int change = start_col_ - start_row_;
    this->template get_arg<0>().set_view(
        bottom_diagonal + change, top_diagonal + change,
        bottom_zero_diagonal + change, top_zero_diagonal + change);
  }

  /**
   * Determine index of bottom diagonal written.
   * @return number of columns
   */
  inline int bottom_diagonal() const {
    return this->template get_arg<0>().bottom_diagonal() - start_col_
           + start_row_;
  }

  /**
   * Determine index of top diagonal written.
   * @return number of columns
   */
  inline int top_diagonal() const {
    return this->template get_arg<0>().top_diagonal() - start_col_ + start_row_;
  }

  /**
   * Evaluates an expression and assigns it to the block.
   * @tparam T_expression type of expression
   * @param rhs input expression
   */
  template <typename T_expression,
            typename
            = require_all_valid_expressions_and_none_scalar_t<T_expression>>
  const block_<T>& operator=(T_expression&& rhs) const {
    auto expression = as_operation_cl(std::forward<T_expression>(rhs));
    if (rows_ * cols_ == 0) {
      return *this;
    }
    expression.evaluate_into(*this);
    return *this;
  }

  /**
   * Checks if desired dimensions match dimensions of the block.
   * @param rows desired number of rows
   * @param cols desired number of columns
   * @throws std::invalid_argument desired dimensions do not match dimensions
   * of the block.
   */
  inline void check_assign_dimensions(int rows, int cols) const {
    check_size_match("block_.check_assign_dimensions", "Rows of ", "block",
                     rows_, "rows of ", "expression", rows);
    check_size_match("block_.check_assign_dimensions", "Columns of ", "block",
                     cols_, "columns of ", "expression", cols);
  }
};

/**
 * Block of a kernel generator expression.
 *
 * Block operation modifies how its argument is indexed. If a matrix is both an
 * argument and result of such an operation (such as in <code> block(a, row1,
 * col1, rows, cols) = block(a, row2, col2, rows, cols);
 * </code>), the result can be wrong due to aliasing. In such case the
 * expression should be evaluating in a temporary by doing <code> block(a, row1,
 * col1, rows, cols) = block(a, row2, col2, rows, cols).eval();</code>. This is
 * not necessary if the bolcks do not overlap or if they are the same block.
 * @tparam T type of argument
 * @param a input argument
 * @param start_row first row of block
 * @param start_col first column of a block
 * @param rows number of rows in block
 * @param cols number of columns in block
 * @return Block of given expression
 */
template <typename T,
          typename = require_all_valid_expressions_and_none_scalar_t<T>>
inline auto block(T&& a, int start_row, int start_col, int rows, int cols) {
  auto&& a_operation = as_operation_cl(std::forward<T>(a)).deep_copy();
  return block_<std::remove_reference_t<decltype(a_operation)>>(
      std::move(a_operation), start_row, start_col, rows, cols);
}

}  // namespace math
}  // namespace stan

#endif
#endif
