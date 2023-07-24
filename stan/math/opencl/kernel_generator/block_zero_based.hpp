#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_BLOCK_ZERO_BASED_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_BLOCK_ZERO_BASED_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl_lhs.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
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
  using base::var_name_;
  using view_transitivity = std::tuple<std::true_type>;
  using base::operator=;

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
    if (start_col < 0) {
      invalid_argument("block", "start_col", start_col,
                       " should be non-negative, but is ");
    }
    if (start_row < 0) {
      invalid_argument("block", "start_row", start_row,
                       " should be non-negative, but is ");
    }
    if (rows < 0) {
      invalid_argument("block", "rows", rows,
                       " should be non-negative, but is ");
    }
    if (cols < 0) {
      invalid_argument("block", "cols", cols,
                       " should be non-negative, but is ");
    }
    if ((a.rows() != base::dynamic && (start_row + rows) > a.rows())
        || (a.cols() != base::dynamic && (start_col + cols) > a.cols())) {
      invalid_argument("block", "block of \"a\"", " is out of bounds", "");
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
   * Generates kernel code for this and nested expressions.
   * @param[in,out] generated map from (pointer to) already generated local
   * operations to variable names
   * @param[in,out] generated_all map from (pointer to) already generated all
   * operations to variable names
   * @param name_gen name generator for this kernel
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether caller already handled matrix view
   * @return part of kernel with code for this and nested expressions
   */
  inline kernel_parts get_kernel_parts(
      std::unordered_map<const void*, const char*>& generated,
      std::unordered_map<const void*, const char*>& generated_all,
      name_generator& name_gen, const std::string& row_index_name,
      const std::string& col_index_name, bool view_handled) const {
    kernel_parts res{};
    if (generated.count(this) == 0) {
      this->var_name_ = name_gen.generate();
      generated[this] = "";
      std::string row_index_name_arg = row_index_name;
      std::string col_index_name_arg = col_index_name;
      modify_argument_indices(row_index_name_arg, col_index_name_arg);
      std::unordered_map<const void*, const char*> generated2;
      res = this->template get_arg<0>().get_kernel_parts(
          generated2, generated_all, name_gen, row_index_name_arg,
          col_index_name_arg, view_handled);
      kernel_parts my
          = this->generate(row_index_name, col_index_name, view_handled,
                           this->template get_arg<0>().var_name_);
      if (generated_all.count(this) == 0) {
        generated_all[this] = "";
      } else {
        my.args = "";
      }
      res += my;
      res.body = res.body_prefix + res.body;
      res.body_prefix = "";
    }
    return res;
  }

  /**
   * Generates kernel code for this expression.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether whether caller already handled matrix view
   * @param var_name_arg name of the variable in kernel that holds argument to
   * this expression
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& row_index_name,
                               const std::string& col_index_name,
                               const bool view_handled,
                               const std::string& var_name_arg) const {
    kernel_parts res;
    res.body
        = type_str<Scalar>() + " " + var_name_ + " = " + var_name_arg + ";\n";
    res.args = "int " + var_name_ + "_i, int " + var_name_ + "_j, ";
    return res;
  }

  /**
   * Sets offset of block to indices of the argument expression
   * @param[in, out] row_index_name row index
   * @param[in, out] col_index_name column index
   */
  inline void modify_argument_indices(std::string& row_index_name,
                                      std::string& col_index_name) const {
    row_index_name = "(" + row_index_name + " + " + var_name_ + "_i)";
    col_index_name = "(" + col_index_name + " + " + var_name_ + "_j)";
  }

  /**
   * Generates kernel code for this expression if it appears on the left hand
   * side of an assignment.
   * @param[in,out] generated map from (pointer to) already generated local
   * operations to variable names
   * @param[in,out] generated_all map from (pointer to) already generated all
   * operations to variable names
   * @param name_gen name generator for this kernel
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @return part of kernel with code for this expressions
   */
  inline kernel_parts get_kernel_parts_lhs(
      std::unordered_map<const void*, const char*>& generated,
      std::unordered_map<const void*, const char*>& generated_all,
      name_generator& name_gen, const std::string& row_index_name,
      const std::string& col_index_name) const {
    if (generated.count(this) == 0) {
      generated[this] = "";
      this->var_name_ = name_gen.generate();
    }
    std::string row_index_name_arg = row_index_name;
    std::string col_index_name_arg = col_index_name;
    modify_argument_indices(row_index_name_arg, col_index_name_arg);
    std::unordered_map<const void*, const char*> generated2;
    kernel_parts res = this->template get_arg<0>().get_kernel_parts_lhs(
        generated2, generated_all, name_gen, row_index_name_arg,
        col_index_name_arg);
    res += this->derived().generate_lhs(row_index_name, col_index_name,
                                        this->template get_arg<0>().var_name_);
    if (generated_all.count(this) == 0) {
      generated_all[this] = "";
    } else {
      res.args = "";
    }
    return res;
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
    res.args = "int " + var_name_ + "_i, int " + var_name_ + "_j, ";
    return res;
  }

  /**
   * Sets kernel arguments for this and nested expressions.
   * @param[in,out] generated map from (pointer to) already generated local
   * operations to variable names
   * @param[in,out] generated_all map from (pointer to) already generated all
   * operations to variable names
   * @param kernel kernel to set arguments on
   * @param[in,out] arg_num consecutive number of the first argument to set.
   * This is incremented for each argument set by this function.
   */
  inline void set_args(
      std::unordered_map<const void*, const char*>& generated,
      std::unordered_map<const void*, const char*>& generated_all,
      cl::Kernel& kernel, int& arg_num) const {
    if (generated.count(this) == 0) {
      generated[this] = "";
      std::unordered_map<const void*, const char*> generated2;
      this->template get_arg<0>().set_args(generated2, generated_all, kernel,
                                           arg_num);
      if (generated_all.count(this) == 0) {
        generated_all[this] = "";
        kernel.setArg(arg_num++, start_row_);
        kernel.setArg(arg_num++, start_col_);
      }
    }
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
    auto& a = this->template get_arg<0>();
    a.set_view(
        bottom_diagonal + change, top_diagonal + change,
        (start_col_ == 0 && start_row_ <= 1 && start_row_ + rows_ == a.rows()
                 && start_col_ + cols_ >= std::min(a.rows() - 1, a.cols())
             ? bottom_zero_diagonal
             : bottom_diagonal)
            + change,
        (start_row_ == 0 && start_col_ <= 1 && start_col_ + cols_ == a.cols()
                 && start_row_ + rows_ >= std::min(a.rows(), a.cols() - 1)
             ? top_zero_diagonal
             : top_diagonal)
            + change);
  }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    std::pair<int, int> arg_diags
        = this->template get_arg<0>().extreme_diagonals();
    return {arg_diags.first - start_col_ + start_row_,
            arg_diags.second - start_col_ + start_row_};
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
 * Block of a kernel generator expression. This version uses zero-based  indices
 * in contrast to the 1-based block(), which is exposed to the Stan language.
 *
 * Block operation modifies how its argument is indexed. If a matrix is both an
 * argument and result of such an operation (such as in <code>
 * block_zero_based(a, row1, col1, rows, cols) = block_zero_based(a, row2, col2,
 * rows, cols);
 * </code>), the result can be wrong due to aliasing. In such case the
 * expression should be evaluating in a temporary by doing <code>
 * block_zero_based(a, row1, col1, rows, cols) = block_zero_based(a, row2, col2,
 * rows, cols).eval();</code>. This is not necessary if the blocks do not
 * overlap or if they are the same block.
 * @tparam T type of argument
 * @param a input argument
 * @param start_row first row of block
 * @param start_col first column of a block
 * @param rows number of rows in block
 * @param cols number of columns in block
 * @return Block of given expression
 */
template <typename T,
          typename = require_all_kernel_expressions_and_none_scalar_t<T>>
inline auto block_zero_based(T&& a, int start_row, int start_col, int rows,
                             int cols) {
  auto&& a_operation = as_operation_cl(std::forward<T>(a)).deep_copy();
  return block_<std::remove_reference_t<decltype(a_operation)>>(
      std::move(a_operation), start_row, start_col, rows, cols);
}
/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
