#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_CONSTANT_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_CONSTANT_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <limits>
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
 * Represents a matrix of single repeated value in kernel generator expressions.
 * @tparam T type of the scalar
 */
template <typename T>
class constant_ : public operation_cl<constant_<T>, T> {
  T a_;
  int rows_;
  int cols_;

 public:
  static_assert(std::is_arithmetic<T>::value,
                "class constant_<T>: std::is_arithmetic<T> must be true!");
  using Scalar = T;
  using base = operation_cl<constant_<T>, T>;
  using base::var_name_;

  /**
   * Constructor
   * @param a scalar value
   * @param rows number of rows of the matrix
   * @param cols number of columns of the matrix
   */
  explicit constant_(const T a, int rows, int cols)
      : a_(a), rows_(rows), cols_(cols) {}

  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline constant_<T> deep_copy() const {
    return constant_<T>(a_, rows_, cols_);
  }

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
    res.args = type_str<Scalar>() + " " + var_name_ + ", ";
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
      kernel.setArg(arg_num++, a_);
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
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    return {std::numeric_limits<int>::min(), std::numeric_limits<int>::max()};
  }
};

/**
 * Matrix of repeated values in kernel generator expressions.
 *
 * In most cases scalars should be directly used instead of this. This is,
 * however, useful for initializing some expression to specific value if that
 * expresssion could also be plain `matrix_cl`.
 *
 * @tparam T type of argument
 * @param a input argument
 * @param rows number of rows
 * @param cols number of columns
 * @return Block of given expression
 */
template <typename T, typename = require_arithmetic_t<T>>
inline auto constant(const T a, int rows, int cols) {
  return constant_<T>(a, rows, cols);
}

/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
