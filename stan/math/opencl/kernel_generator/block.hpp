#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_BLOCK_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_BLOCK_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/err/throw_domain_error.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl_lhs.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/is_valid_expression.hpp>
#include <string>
#include <type_traits>
#include <set>
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
  using base::arguments_;

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
   * generates kernel code for this and nested expressions.
   * @param[in,out] generated set of already generated operation_cls
   * @param ng name generator for this kernel
   * @param i row index variable name
   * @param j column index variable name
   * @return part of kernel with code for this and nested expressions
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
   * generates kernel code for this and nested expressions if this expression
   * appears on the left hand side of an assignment.
   * @param[in,out] generated set of already generated operation_cls
   * @param ng name generator for this kernel
   * @param i row index variable name
   * @param j column index variable name
   * @return part of kernel with code for this and nested expressions
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
      std::get<0>(arguments_).set_args(generated, kernel, arg_num);
      kernel.setArg(arg_num++, start_row_);
      kernel.setArg(arg_num++, start_col_);
    }
  }

  /**
   * View of a matrix that would be the result of evaluating this expression.
   * @return view
   */
  inline matrix_cl_view view() const { return std::get<0>(arguments_).view(); }

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
   * Evaluates an expression and assigns it to the block.
   * @tparam T_expression type of expression
   * @param rhs input expression
   */
  template <typename T_expression,
            typename
            = require_all_valid_expressions_and_none_scalar_t<T_expression>>
  const block_<T>& operator=(T_expression&& rhs) const {
    check_size_match("block.operator=", "Rows of ", "rhs", rhs.rows(),
                     "rows of ", "*this", this->rows());
    check_size_match("block.operator=", "Cols of ", "rhs", rhs.cols(),
                     "cols of ", "*this", this->cols());
    auto expression = as_operation_cl(std::forward<T_expression>(rhs));
    expression.evaluate_into(*this);
    return *this;
  }
};

/**
 * Block of a kernel generator expression.
 * @tparam T type of argument
 * @param a input argument
 * @return Block of given expression
 */
template <typename T,
          typename = require_all_valid_expressions_and_none_scalar_t<T>>
inline block_<as_operation_cl_t<T>> block(T&& a, int start_row, int start_col,
                                          int rows, int cols) {
  return block_<as_operation_cl_t<T>>(as_operation_cl(std::forward<T>(a)),
                                      start_row, start_col, rows, cols);
}

}  // namespace math
}  // namespace stan

#endif
#endif
