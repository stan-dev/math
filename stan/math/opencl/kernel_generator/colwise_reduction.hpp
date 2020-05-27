#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_COLWISE_REDUCTION_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_COLWISE_REDUCTION_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/is_kernel_expression.hpp>
#include <stan/math/opencl/kernel_generator/rowwise_reduction.hpp>
#include <set>
#include <string>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {
/** \addtogroup opencl_kernel_generator
 *  @{
 */

/**
 * Represents a column wise reduction in kernel generator expressions. So as to
 * be efficient column wise reductions are only done partially. That means
 * instead of 1 row kernel output will have a few rows that need to be reduced
 * to obtain final result. This can be done in a separate kernel or after
 * copying to CPU. Also column wise reductions can not be used as arguments to
 * other operations - they can only be evaluated.
 * @tparam Derived derived type
 * @tparam T type of first argument
 * @tparam operation type with member function generate that accepts two
 * variable names and returns OpenCL source code for reduction operation_cl
 * @tparam PassZero whether \c operation passes trough zeros
 */
template <typename Derived, typename T, typename Operation>
class colwise_reduction
    : public operation_cl<Derived, typename std::remove_reference_t<T>::Scalar,
                          T> {
 public:
  using Scalar = typename std::remove_reference_t<T>::Scalar;
  using base = operation_cl<Derived, Scalar, T>;
  using base::var_name_;
  static const bool require_specific_local_size = true;

 protected:
  std::string init_;
  using base::derived;

 public:
  using base::cols;
  /**
   * Constructor
   * @param a the expression to reduce
   * @param init OpenCL source code of initialization value for reduction
   */
  explicit colwise_reduction(T&& a, const std::string& init)
      : base(std::forward<T>(a)), init_(init) {}

  /**
   * Generates kernel code for assigning this expression into result expression.
   * @param[in,out] generated set of (pointer to) already generated operations
   * @param ng name generator for this kernel
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param result expression into which result is to be assigned
   * @return part of kernel with code for this and nested expressions
   */
  template <typename T_result>
  kernel_parts get_whole_kernel_parts(
      std::set<const operation_cl_base*>& generated, name_generator& ng,
      const std::string& row_index_name, const std::string& col_index_name,
      const T_result& result) const {
    kernel_parts parts = derived().get_kernel_parts(
        generated, ng, row_index_name, col_index_name, false);
    kernel_parts out_parts = result.get_kernel_parts_lhs(
        generated, ng, row_index_name, col_index_name);

    parts.args += out_parts.args;
    parts.reduction += "if (lid_i == 0) {\n"
                     + result.var_name_
                     + "_global[j * blocks_rows + wg_id_i] = "
                     + derived().var_name_ + "_local[0];\n"
                     "}\n";
    return parts;
  }

  /**
   * Generates kernel code for this and nested expressions.
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param view_handled whether whether caller already handled matrix view
   * @param var_name_arg name of the variable in kernel that holds argument to
   * this expression
   * @return part of kernel with code for this and nested expressions
   */
  inline kernel_parts generate(const std::string& row_index_name,
                               const std::string& col_index_name,
                               const bool view_handled,
                               const std::string& var_name_arg) const {
    kernel_parts res;
    res.initialization = type_str<Scalar>() + " " + var_name_ + " = " + init_
        + ";\n" "__local " + type_str<Scalar>() + " "
        + var_name_ + "_local[LOCAL_SIZE_];\n";
    res.body = var_name_ + " = " + var_name_arg + ";\n";
    res.reduction =
          var_name_ + "_local[lid_i] = " + var_name_ + ";\n"
          "barrier(CLK_LOCAL_MEM_FENCE);\n"
          "for (int step = lsize_i / REDUCTION_STEP_SIZE; "
                "step > 0; step /= REDUCTION_STEP_SIZE) {\n"
          "  if (lid_i < step) {\n"
          "    for (int i = 1; i < REDUCTION_STEP_SIZE; i++) {\n"
          "      " + var_name_ + "_local[lid_i] = " +
        Operation::generate(var_name_ + "_local[lid_i]",
                            var_name_ + "_local[lid_i + step * i]") + ";\n"
          "    }\n"
          "  }\n"
          "  barrier(CLK_LOCAL_MEM_FENCE);\n"
          "}\n";
    return res;
  }

  /**
   * Number of rows of a matrix that would be the result of evaluating this
   * expression.
   * @return number of rows
   */
  inline int rows() const {
    int local_rows = opencl_context.base_opts().at("LOCAL_SIZE_");
    int wgs_rows
        = (this->template get_arg<0>().rows() + local_rows - 1) / local_rows;
    return wgs_rows;
  }

  /**
   * Number of rows threads need to be launched for.
   * @return number of rows
   */
  inline int thread_rows() const { return this->template get_arg<0>().rows(); }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    return {-rows() + 1, cols() - 1};
  }
};  // namespace math

/**
 * Represents column wise sum - reduction in kernel generator expressions.
 * @tparam T type of expression
 */
template <typename T>
class colwise_sum_ : public colwise_reduction<colwise_sum_<T>, T, sum_op> {
  using base = colwise_reduction<colwise_sum_<T>, T, sum_op>;
  using base::arguments_;

 public:
  explicit colwise_sum_(T&& a)
      : colwise_reduction<colwise_sum_<T>, T, sum_op>(std::forward<T>(a), "0") {
  }
  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return colwise_sum_<std::remove_reference_t<decltype(arg_copy)>>(
        std::move(arg_copy));
  }
};

/**
 * Column wise sum - reduction of a kernel generator expression. So as to
 * be efficient column wise reductions are only done partially. That means
 * instead of 1 row kernel output will have a few rows that need to be reduced
 * to obtain final result. This can be done in a separate kernel or after
 * copying to CPU. Also column wise reductions can not be used as arguments to
 * other operations - they can only be evaluated.
 * @tparam T type of input expression
 * @param a expression to reduce
 * @return sum
 */
template <typename T,
          typename = require_all_kernel_expressions_and_none_scalar_t<T>>
inline auto colwise_sum(T&& a) {
  auto&& arg_copy = as_operation_cl(std::forward<T>(a)).deep_copy();
  return colwise_sum_<std::remove_reference_t<decltype(arg_copy)>>(
      std::move(arg_copy));
}

/**
 * Represents column wise max - reduction in kernel generator expressions.
 * @tparam T type of expression
 */
template <typename T>
class colwise_max_ : public colwise_reduction<
                         colwise_max_<T>, T,
                         max_op<typename std::remove_reference_t<T>::Scalar>> {
  using base
      = colwise_reduction<colwise_max_<T>, T,
                          max_op<typename std::remove_reference_t<T>::Scalar>>;
  using base::arguments_;

 public:
  using op = max_op<typename std::remove_reference_t<T>::Scalar>;
  explicit colwise_max_(T&& a)
      : colwise_reduction<colwise_max_<T>, T, op>(std::forward<T>(a),
                                                  op::init()) {}
  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return colwise_max_<std::remove_reference_t<decltype(arg_copy)>>(
        std::move(arg_copy));
  }
};

/**
 * Column wise max - reduction of a kernel generator expression. So as to
 * be efficient column wise reductions are only done partially. That means
 * instead of 1 row kernel output will have a few rows that need to be reduced
 * to obtain final result. This can be done in a separate kernel or after
 * copying to CPU. Also column wise reductions can not be used as arguments to
 * other operations - they can only be evaluated.
 * @tparam T type of input expression
 * @param a expression to reduce
 * @return max
 */
template <typename T,
          typename = require_all_kernel_expressions_and_none_scalar_t<T>>
inline auto colwise_max(T&& a) {
  auto&& arg_copy = as_operation_cl(std::forward<T>(a)).deep_copy();
  return colwise_max_<std::remove_reference_t<decltype(arg_copy)>>(
      std::move(arg_copy));
}

/**
 * Represents column wise min - reduction in kernel generator expressions.
 * @tparam T type of expression
 */
template <typename T>
class colwise_min_ : public colwise_reduction<
                         colwise_min_<T>, T,
                         min_op<typename std::remove_reference_t<T>::Scalar>> {
  using base
      = colwise_reduction<colwise_min_<T>, T,
                          min_op<typename std::remove_reference_t<T>::Scalar>>;
  using base::arguments_;

 public:
  using op = min_op<typename std::remove_reference_t<T>::Scalar>;
  explicit colwise_min_(T&& a)
      : colwise_reduction<colwise_min_<T>, T, op>(std::forward<T>(a),
                                                  op::init()) {}
  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return colwise_min_<std::remove_reference_t<decltype(arg_copy)>>(
        std::move(arg_copy));
  }
};

/**
 * Column wise min - reduction of a kernel generator expression.  So as to
 * be efficient column wise reductions are only done partially. That means
 * instead of 1 row kernel output will have a few rows that need to be reduced
 * to obtain final result. This can be done in a separate kernel or after
 * copying to CPU. Also column wise reductions can not be used as arguments to
 * other operations - they can only be evaluated.
 * @tparam T type of input expression
 * @param a expression to reduce
 * @return min
 */
template <typename T,
          typename = require_all_kernel_expressions_and_none_scalar_t<T>>
inline auto colwise_min(T&& a) {
  auto&& arg_copy = as_operation_cl(std::forward<T>(a)).deep_copy();
  return colwise_min_<std::remove_reference_t<decltype(arg_copy)>>(
      std::move(arg_copy));
}
/** @}*/
}  // namespace math
}  // namespace stan
#endif
#endif
