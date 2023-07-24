#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_APPEND_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_APPEND_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/err/check_opencl.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/scalar.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/common_return_scalar.hpp>
#include <algorithm>
#include <string>
#include <tuple>
#include <type_traits>
#include <map>
#include <utility>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */

/**
 * Represents appending of rows in kernel generator expressions.
 * @tparam T_a type of first argument
 * @tparam T_b type of second argument
 */
template <typename T_a, typename T_b>
class append_row_ : public operation_cl<append_row_<T_a, T_b>,
                                        common_scalar_t<T_a, T_b>, T_a, T_b> {
 public:
  using Scalar = common_scalar_t<T_a, T_b>;
  using base = operation_cl<append_row_<T_a, T_b>, Scalar, T_a, T_b>;
  using base::var_name_;

 protected:
  using base::arguments_;

 public:
  /**
   * Constructor
   * @param a first argument
   * @param b second argument
   */
  append_row_(T_a&& a, T_b&& b)  // NOLINT
      : base(std::forward<T_a>(a), std::forward<T_b>(b)) {
    if (a.cols() != base::dynamic && b.cols() != base::dynamic) {
      check_size_match("append_row", "Columns of ", "a", a.cols(),
                       "columns of ", "b", b.cols());
    }
    if (a.rows() < 0) {
      invalid_argument("append_row", "Rows of a", a.rows(),
                       "should be nonnegative!");
    }
    if (b.rows() < 0) {
      invalid_argument("append_row", "Rows of b", b.rows(),
                       "should be nonnegative!");
    }
  }

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& a_copy = this->template get_arg<0>().deep_copy();
    auto&& b_copy = this->template get_arg<1>().deep_copy();
    return append_row_<std::remove_reference_t<decltype(a_copy)>,
                       std::remove_reference_t<decltype(b_copy)>>{
        std::move(a_copy), std::move(b_copy)};
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
      var_name_ = name_gen.generate();
      generated[this] = "";
      kernel_parts parts_a = this->template get_arg<0>().get_kernel_parts(
          generated, generated_all, name_gen, row_index_name, col_index_name,
          true);
      std::string row_index_name_b
          = "(" + row_index_name + " - " + var_name_ + "_first_rows)";
      std::unordered_map<const void*, const char*> generated_b;
      kernel_parts parts_b = this->template get_arg<1>().get_kernel_parts(
          generated_b, generated_all, name_gen, row_index_name_b,
          col_index_name, true);
      res = parts_a + parts_b;
      res.body = type_str<Scalar>() + " " + var_name_ + ";\n"
          "if("+ row_index_name +" < " + var_name_ + "_first_rows){\n"
          + parts_a.body +
          var_name_ + " = " + this->template get_arg<0>().var_name_ + ";\n"
          "} else{\n"
          + parts_b.body +
          var_name_ + " = " + this->template get_arg<1>().var_name_ + ";\n"
          "}\n";
      res.args += "int " + var_name_ + "_first_rows, ";
    }
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
      this->template get_arg<0>().set_args(generated, generated_all, kernel,
                                           arg_num);
      std::unordered_map<const void*, const char*> generated_b;
      this->template get_arg<1>().set_args(generated_b, generated_all, kernel,
                                           arg_num);
      kernel.setArg(arg_num++, this->template get_arg<0>().rows());
    }
  }

  /**
   * Number of rows of a matrix that would be the result of evaluating this
   * expression.
   * @return number of rows
   */
  inline int rows() const {
    return this->template get_arg<0>().rows()
           + this->template get_arg<1>().rows();
  }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    std::pair<int, int> a_diags
        = this->template get_arg<0>().extreme_diagonals();
    std::pair<int, int> b_diags
        = this->template get_arg<1>().extreme_diagonals();
    int my_rows = this->template get_arg<0>().rows();
    return {std::min(a_diags.first, b_diags.first - my_rows),
            std::max(a_diags.second, b_diags.second - my_rows)};
  }
};

/**
 * Stack the rows of the first argument on top of the second argument.
 *
 * @tparam Ta type of first argument
 * @tparam Ta type of second argument
 * @param a First argument
 * @param b Second argument
 * @return Stacked arguments
 */
template <typename Ta, typename Tb,
          typename = require_all_kernel_expressions_and_none_scalar_t<Ta, Tb>>
inline auto append_row(Ta&& a, Tb&& b) {
  auto&& a_operation = as_operation_cl(std::forward<Ta>(a)).deep_copy();
  auto&& b_operation = as_operation_cl(std::forward<Tb>(b)).deep_copy();
  return append_row_<std::remove_reference_t<decltype(a_operation)>,
                     std::remove_reference_t<decltype(b_operation)>>(
      std::move(a_operation), std::move(b_operation));
}

/**
 * Represents appending of cols in kernel generator expressions.
 * @tparam T_a type of first argument
 * @tparam T_b type of second argument
 */
template <typename T_a, typename T_b>
class append_col_ : public operation_cl<append_col_<T_a, T_b>,
                                        common_scalar_t<T_a, T_b>, T_a, T_b> {
 public:
  using Scalar = common_scalar_t<T_a, T_b>;
  using base = operation_cl<append_col_<T_a, T_b>, Scalar, T_a, T_b>;
  using base::var_name_;

 protected:
  using base::arguments_;

 public:
  /**
   * Constructor
   * @param a first argument
   * @param b second argument
   */
  append_col_(T_a&& a, T_b&& b)  // NOLINT
      : base(std::forward<T_a>(a), std::forward<T_b>(b)) {
    if (a.rows() != base::dynamic && b.rows() != base::dynamic) {
      check_size_match("append_col", "Rows of ", "a", a.rows(), "rows of ", "b",
                       b.rows());
    }
    if (a.cols() < 0) {
      invalid_argument("append_col", "Columns of a", a.cols(),
                       "should be nonnegative!");
    }
    if (b.cols() < 0) {
      invalid_argument("append_col", "Columns of b", b.cols(),
                       "should be nonnegative!");
    }
  }

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& a_copy = this->template get_arg<0>().deep_copy();
    auto&& b_copy = this->template get_arg<1>().deep_copy();
    return append_col_<std::remove_reference_t<decltype(a_copy)>,
                       std::remove_reference_t<decltype(b_copy)>>{
        std::move(a_copy), std::move(b_copy)};
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
      var_name_ = name_gen.generate();
      generated[this] = "";
      kernel_parts parts_a = this->template get_arg<0>().get_kernel_parts(
          generated, generated_all, name_gen, row_index_name, col_index_name,
          true);
      std::string col_index_name_b
          = "(" + col_index_name + " - " + var_name_ + "_first_cols)";
      std::unordered_map<const void*, const char*> generated_b;
      kernel_parts parts_b = this->template get_arg<1>().get_kernel_parts(
          generated_b, generated_all, name_gen, row_index_name,
          col_index_name_b, true);
      res = parts_a + parts_b;
      res.body = type_str<Scalar>() + " " + var_name_ + ";\n"
          "if("+ col_index_name +" < " + var_name_ + "_first_cols){\n"
          + parts_a.body +
          var_name_ + " = " + this->template get_arg<0>().var_name_ + ";\n"
          "} else{\n"
          + parts_b.body +
          var_name_ + " = " + this->template get_arg<1>().var_name_ + ";\n"
          "}\n";
      res.args += "int " + var_name_ + "_first_cols, ";
    }
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
      this->template get_arg<0>().set_args(generated, generated_all, kernel,
                                           arg_num);
      std::unordered_map<const void*, const char*> generated_b;
      this->template get_arg<1>().set_args(generated_b, generated_all, kernel,
                                           arg_num);
      kernel.setArg(arg_num++, this->template get_arg<0>().cols());
    }
  }

  /**
   * Number of rows of a matrix that would be the result of evaluating this
   * expression.
   * @return number of rows
   */
  inline int cols() const {
    return this->template get_arg<0>().cols()
           + this->template get_arg<1>().cols();
  }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    std::pair<int, int> a_diags
        = this->template get_arg<0>().extreme_diagonals();
    std::pair<int, int> b_diags
        = this->template get_arg<1>().extreme_diagonals();
    int my_cols = this->template get_arg<0>().cols();
    return {std::min(a_diags.first, b_diags.first + my_cols),
            std::max(a_diags.second, b_diags.second + my_cols)};
  }
};

/**
 * Stack the cols of the arguments.
 *
 * @tparam Ta type of first argument
 * @tparam Ta type of second argument
 * @param a First argument
 * @param b Second argument
 * @return Stacked arguments
 */
template <typename Ta, typename Tb,
          typename = require_all_kernel_expressions_and_none_scalar_t<Ta, Tb>>
inline auto append_col(Ta&& a, Tb&& b) {
  auto&& a_operation = as_operation_cl(std::forward<Ta>(a)).deep_copy();
  auto&& b_operation = as_operation_cl(std::forward<Tb>(b)).deep_copy();
  return append_col_<std::remove_reference_t<decltype(a_operation)>,
                     std::remove_reference_t<decltype(b_operation)>>(
      std::move(a_operation), std::move(b_operation));
}

}  // namespace math
}  // namespace stan

#endif
#endif
