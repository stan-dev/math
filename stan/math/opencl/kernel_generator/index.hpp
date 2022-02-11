#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_INDEX_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_INDEX_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/operation_cl.hpp>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */

/**
 * Represents operation that determines row index.
 */
class row_index : public operation_cl<row_index, int> {
 public:
  using Scalar = int;
  using base = operation_cl<row_index, int>;
  using base::var_name_;
  int rows_;
  int cols_;

  /**
   * Constructor for given number of rows and columns.
   * @param rows number of rows of this expression (default = broadcast)
   * @param cols number of columns of this expression (default = broadcast)
   */
  explicit row_index(int rows = base::dynamic, int cols = base::dynamic)
      : rows_(rows), cols_(cols) {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline row_index deep_copy() const { return row_index(rows_, cols_); }

  /**
   * Generates kernel code for this expression.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether whether caller already handled matrix view
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& row_index_name,
                               const std::string& col_index_name,
                               const bool view_handled) const {
    kernel_parts res{};
    res.body = "int " + var_name_ + " = " + row_index_name + ";\n";
    return res;
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
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    return {std::numeric_limits<int>::min(), std::numeric_limits<int>::max()};
  }
};

/**
 * Represents operation that determines column index.
 */
class col_index : public operation_cl<col_index, int> {
 public:
  using Scalar = int;
  using base = operation_cl<col_index, int>;
  using base::var_name_;
  int rows_;
  int cols_;

  /**
   * Constructor for given number of rows and columns.
   * @param rows number of rows of this expression (default = broadcast)
   * @param cols number of columns of this expression (default = broadcast)
   */
  explicit col_index(int rows = base::dynamic, int cols = base::dynamic)
      : rows_(rows), cols_(cols) {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline col_index deep_copy() const { return col_index(rows_, cols_); }

  /**
   * Generates kernel code for this expression.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether whether caller already handled matrix view
   * @return part of kernel with code for this expression
   */
  inline kernel_parts generate(const std::string& row_index_name,
                               const std::string& col_index_name,
                               const bool view_handled) const {
    kernel_parts res{};
    res.body = "int " + var_name_ + " = " + col_index_name + ";\n";
    return res;
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
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    return {std::numeric_limits<int>::min(), std::numeric_limits<int>::max()};
  }
};
/** @}*/

}  // namespace math
}  // namespace stan

#endif
#endif
