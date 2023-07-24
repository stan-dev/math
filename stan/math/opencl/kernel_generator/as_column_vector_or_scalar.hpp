#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_AS_COLUMN_VECTOR_OR_SCALAR_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_AS_COLUMN_VECTOR_OR_SCALAR_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl_lhs.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/constant.hpp>
#include <stan/math/prim/err/check_vector.hpp>
#include <stan/math/prim/meta.hpp>
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
 * Represents as_column_vector_or_scalar of a row or column vector in
 * kernel generator expressions.
 * @tparam Derived derived type
 * @tparam T type of argument
 */
template <typename T>
class as_column_vector_or_scalar_
    : public operation_cl_lhs<as_column_vector_or_scalar_<T>,
                              typename std::remove_reference_t<T>::Scalar, T> {
 public:
  using Scalar = typename std::remove_reference_t<T>::Scalar;
  using base = operation_cl_lhs<as_column_vector_or_scalar_<T>, Scalar, T>;
  using base::var_name_;
  using base::operator=;

  /**
   * Constructor
   * @param a expression (must be row or column vector)
   */
  explicit as_column_vector_or_scalar_(T&& a) : base(std::forward<T>(a)) {
    check_vector("as_column_vector_or_scalar", "a", a);
  }

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return as_column_vector_or_scalar_<
        std::remove_reference_t<decltype(arg_copy)>>{std::move(arg_copy)};
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
    res.args = "int " + var_name_ + "_transpose, ";
    res.body
        = type_str<Scalar>() + " " + var_name_ + " = " + var_name_arg + ";\n";
    return res;
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
    res.args = "int " + var_name_ + "_transpose, ";
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
        kernel.setArg(arg_num++,
                      static_cast<int>(this->template get_arg<0>().rows()
                                       < this->template get_arg<0>().cols()));
      }
    }
  }

  /**
   * Swaps indices \c row_index_name and \c col_index_name for the argument
   * expression if necessary.
   * @param[in, out] row_index_name row index
   * @param[in, out] col_index_name column index
   */
  inline void modify_argument_indices(std::string& row_index_name,
                                      std::string& col_index_name) const {
    std::string row_index_name2 = "(" + var_name_ + "_transpose ? "
                                  + col_index_name + " : " + row_index_name
                                  + ")";
    col_index_name = "(" + var_name_ + "_transpose ? " + row_index_name + " : "
                     + col_index_name + ")";
    row_index_name = std::move(row_index_name2);
  }

  /**
   * Number of rows of a matrix that would be the result of evaluating this
   * expression.
   * @return number of rows
   */
  inline int rows() const {
    return std::max(this->template get_arg<0>().rows(),
                    this->template get_arg<0>().cols());
  }

  /**
   * Number of columns of a matrix that would be the result of evaluating this
   * expression.
   * @return 1
   */
  inline int cols() const {
    return std::min(this->template get_arg<0>().rows(),
                    this->template get_arg<0>().cols());
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
  inline void set_view(int bottom_diagonal, int top_diagonal,
                       int bottom_zero_diagonal, int top_zero_diagonal) const {
    auto& arg = this->template get_arg<0>();
    if (arg.rows() >= arg.cols()) {
      arg.set_view(bottom_diagonal, top_diagonal, bottom_zero_diagonal,
                   top_zero_diagonal);
    } else {
      arg.set_view(top_diagonal, bottom_diagonal, top_zero_diagonal,
                   bottom_zero_diagonal);
    }
  }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    auto& arg = this->template get_arg<0>();
    std::pair<int, int> arg_diags = arg.extreme_diagonals();
    if (arg.rows() >= arg.cols()) {
      return arg_diags;
    } else {
      return {-arg_diags.second, -arg_diags.first};
    }
  }

  /**
   * Checks if desired dimensions match dimensions of the argument vector.
   * @param rows desired number of rows
   * @param cols desired number of columns
   * @throws std::invalid_argument desired dimensions do not match dimensions
   * of the block.
   */
  inline void check_assign_dimensions(int rows, int cols) const {
    // use a dummy expression with same number of rows and cols to simplify the
    // check
    check_vector("as_column_vector_or_scalar_::check_assign_dimensions()",
                 "expression assigned to as_column_vector_or_scalar",
                 constant(0, rows, cols));
    auto& arg = this->template get_arg<0>();
    int arg_rows = arg.rows();
    int arg_cols = arg.cols();
    if (arg_rows >= arg_cols) {
      check_size_match("block_.check_assign_dimensions", "Rows of ",
                       "check_assign_dimensions argument", arg_rows, "rows of ",
                       "expression", rows);
      check_size_match("block_.check_assign_dimensions", "Columns of ",
                       "check_assign_dimensions argument", arg_cols,
                       "columns of ", "expression", cols);
    } else {
      check_size_match("block_.check_assign_dimensions", "Columns of ",
                       "check_assign_dimensions argument", arg_cols, "rows of ",
                       "expression", rows);
      check_size_match("block_.check_assign_dimensions", "Rows of ",
                       "check_assign_dimensions argument", arg_rows,
                       "columns of ", "expression", cols);
    }
  }
};

/**
 * as_column_vector_or_scalar of a kernel generator expression.
 *
 * @tparam T type of argument
 * @param a input argument (must be a row or a column vector)
 * @return as_column_vector_or_scalar of given expression
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline auto as_column_vector_or_scalar(T&& a) {
  auto&& a_operation = as_operation_cl(std::forward<T>(a)).deep_copy();
  return as_column_vector_or_scalar_<
      std::remove_reference_t<decltype(a_operation)>>(std::move(a_operation));
}
/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
