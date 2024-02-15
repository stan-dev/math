#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_REDUCTION_2D_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_REDUCTION_2D_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/opencl_context.hpp>
#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/colwise_reduction.hpp>
#include <stan/math/opencl/kernel_generator/rowwise_reduction.hpp>
#include <stan/math/opencl/kernel_generator/calc_if.hpp>
#include <map>
#include <string>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {
/** \addtogroup opencl_kernel_generator
 *  @{
 */

namespace internal {
class reduction_2d_base {};
}  // namespace internal

/**
 * Represents a two dimensional reduction in kernel generator expressions. So as
 * to be efficient two dimensional reductions are only done partially. That
 * means instead of 1 element kernel output can have a few rows and a few
 * columns that need to be reduced to obtain final result (actually it is 1
 * result per work group run - roughly 16 times the number of compute units on
 * the OpenCL device). This can be done in a separate kernel or after copying to
 * CPU. Also two dimensional reductions can not be used as arguments to other
 * operations - they can only be evaluated.
 * @tparam Derived Type derived from `reduction_2d`
 * @tparam T type of first argument
 * @tparam Operation type with member function generate that accepts two
 * variable names and returns OpenCL source code for reduction operation_cl
 */
template <typename Derived, typename T, typename Operation>
class reduction_2d
    : public internal::reduction_2d_base,
      public operation_cl<Derived, typename std::remove_reference_t<T>::Scalar,
                          T> {
 public:
  using Scalar = typename std::remove_reference_t<T>::Scalar;
  using base = operation_cl<Derived, Scalar, T>;
  using base::var_name_;
  static constexpr bool require_specific_local_size = true;

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
  explicit reduction_2d(T&& a, const std::string& init)
      : base(std::forward<T>(a)), init_(init) {}

  /**
   * Generates kernel code for assigning this expression into result expression.
   * @param[in,out] generated map from (pointer to) already generated local
   * operations to variable names
   * @param[in,out] generated_all map from (pointer to) already generated all
   * operations to variable names
   * @param ng name generator for this kernel
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param result expression into which result is to be assigned
   * @return part of kernel with code for this and nested expressions
   */
  template <typename T_result>
  kernel_parts get_whole_kernel_parts(
      std::unordered_map<const void*, const char*>& generated,
      std::unordered_map<const void*, const char*>& generated_all,
      name_generator& ng, const std::string& row_index_name,
      const std::string& col_index_name, const T_result& result) const {
    kernel_parts parts = derived().get_kernel_parts(
        generated, generated_all, ng, row_index_name, col_index_name, false);
    kernel_parts out_parts = result.get_kernel_parts_lhs(
        generated, generated_all, ng, row_index_name, col_index_name);

    parts.args += out_parts.args;
    parts.reduction_2d += "if (lid_i == 0) {\n"
                     + result.var_name_
                     + "_global[wg_id_j * n_groups_i + wg_id_i] = "
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
    res.declarations = "__local " + type_str<Scalar>() + " " + var_name_
                       + "_local[LOCAL_SIZE_];\n" + type_str<Scalar>() + " "
                       + var_name_ + " = " + init_ + ";\n";
    res.body = var_name_ + " = " + Operation::generate(var_name_, var_name_arg)
               + ";\n";
    res.reduction_2d =
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
    int arg_rows = this->template get_arg<0>().rows();
    int arg_cols = this->template get_arg<0>().cols();
    if (arg_cols == 0) {
      return 1;
    }
    if (arg_cols == base::dynamic) {
      return base::dynamic;
    }
    return internal::colwise_reduction_wgs_rows(arg_rows, arg_cols);
  }

  /**
   * Number of columns of a matrix that would be the result of evaluating this
   * expression.
   * @return number of columns
   */
  inline int cols() const {
    int arg_rows = this->template get_arg<0>().rows();
    int arg_cols = this->template get_arg<0>().cols();
    if (arg_cols == 0) {
      return 0;
    }
    if (arg_cols == base::dynamic) {
      return base::dynamic;
    }
    int wgs_rows = internal::colwise_reduction_wgs_rows(arg_rows, arg_cols);
    if (wgs_rows == 0) {
      return 0;
    }
    return (arg_cols + wgs_rows - 1) / wgs_rows;
  }

  /**
   * Number of rows threads need to be launched for.
   * @return number of rows
   */
  inline int thread_rows() const { return this->template get_arg<0>().rows(); }

  /**
   * Number of rows threads need to be launched for.
   * @return number of rows
   */
  inline int thread_cols() const { return this->template get_arg<0>().cols(); }

  /**
   * Determine indices of extreme sub- and superdiagonals written.
   * @return pair of indices - bottom and top diagonal
   */
  inline std::pair<int, int> extreme_diagonals() const {
    return {-rows() + 1, cols() - 1};
  }
};  // namespace math

/**
 * Represents two dimensional sum - reduction in kernel generator expressions.
 * @tparam T type of expression
 */
template <typename T>
class sum_2d_ : public reduction_2d<sum_2d_<T>, T, sum_op> {
  using base = reduction_2d<sum_2d_<T>, T, sum_op>;
  using base::arguments_;

 public:
  explicit sum_2d_(T&& a)
      : reduction_2d<sum_2d_<T>, T, sum_op>(std::forward<T>(a), "0") {}
  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return sum_2d_<std::remove_reference_t<decltype(arg_copy)>>(
        std::move(arg_copy));
  }
};

/**
 * Two dimensional sum - reduction of a kernel generator expression. So as to
 * be efficient two dimensional reductions are only done partially. That means
 * instead of 1 element kernel output can have a few rows and a few columns
 * that need to be reduced to obtain final result (actually it is 1 result per
 * work group run - roughly 16 times the number of compute units on the OpenCL
 * device). This can be done in a separate kernel or after copying to CPU. Also
 * two dimensional reductions can not be used as arguments to other operations -
 * they can only be evaluated.
 * @tparam T type of input expression
 * @param a the expression to reduce
 * @return sum
 */
template <typename T, require_all_kernel_expressions_t<T>* = nullptr>
inline auto sum_2d(T&& a) {
  auto&& arg_copy = as_operation_cl(std::forward<T>(a)).deep_copy();
  return sum_2d_<as_operation_cl_t<T>>(as_operation_cl(std::forward<T>(a)));
}

/**
 * Represents two dimensional product - reduction in kernel generator
 * expressions.
 * @tparam T type of expression
 */
template <typename T>
class prod_2d_ : public reduction_2d<prod_2d_<T>, T, prod_op> {
  using base = reduction_2d<prod_2d_<T>, T, prod_op>;
  using base::arguments_;

 public:
  explicit prod_2d_(T&& a)
      : reduction_2d<prod_2d_<T>, T, prod_op>(std::forward<T>(a), "1") {}
  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return prod_2d_<std::remove_reference_t<decltype(arg_copy)>>(
        std::move(arg_copy));
  }
};

/**
 * Two dimensional product - reduction of a kernel generator expression. So as
 * to be efficient two dimensional reductions are only done partially. That
 * means instead of 1 element kernel output can have a few rows and a few
 * columns that need to be reduced to obtain final result (actually it is 1
 * result per work group run - roughly 16 times the number of compute units on
 * the OpenCL device). This can be done in a separate kernel or after copying to
 * CPU. Also two dimensional reductions can not be used as arguments to other
 * operations - they can only be evaluated.
 * @tparam T type of input expression
 * @param a the expression to reduce
 * @return prod
 */
template <typename T, require_all_kernel_expressions_t<T>* = nullptr>
inline auto prod_2d(T&& a) {
  auto&& arg_copy = as_operation_cl(std::forward<T>(a)).deep_copy();
  return prod_2d_<as_operation_cl_t<T>>(as_operation_cl(std::forward<T>(a)));
}

/**
 * Represents two dimensional max - reduction in kernel generator expressions.
 * @tparam T type of expression
 */
template <typename T>
class max_2d_
    : public reduction_2d<max_2d_<T>, T,
                          max_op<typename std::remove_reference_t<T>::Scalar>> {
  using base
      = reduction_2d<max_2d_<T>, T,
                     max_op<typename std::remove_reference_t<T>::Scalar>>;
  using base::arguments_;

 public:
  using op = max_op<typename std::remove_reference_t<T>::Scalar>;
  explicit max_2d_(T&& a)
      : reduction_2d<max_2d_<T>, T, op>(std::forward<T>(a), op::init()) {}
  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return max_2d_<std::remove_reference_t<decltype(arg_copy)>>(
        std::move(arg_copy));
  }
};

/**
 * Two dimensional max - reduction of a kernel generator expression. So as to
 * be efficient two dimensional reductions are only done partially. That means
 * instead of 1 element kernel output can have a few rows and a few columns that
 * need to be reduced to obtain final result (actually it is 1 result per work
 * group run - roughly 16 times the number of compute units on the OpenCL
 * device). This can be done in a separate kernel or after copying to CPU. Also
 * two dimensional reductions can not be used as arguments to other operations -
 * they can only be evaluated.
 * @tparam T type of input expression
 * @param a the expression to reduce
 * @return max
 */
template <typename T, require_all_kernel_expressions_t<T>* = nullptr>
inline auto max_2d(T&& a) {
  auto&& arg_copy = as_operation_cl(std::forward<T>(a)).deep_copy();
  return max_2d_<as_operation_cl_t<T>>(as_operation_cl(std::forward<T>(a)));
}

/**
 * Represents two dimensional min - reduction in kernel generator expressions.
 * @tparam T type of expression
 */
template <typename T>
class min_2d_
    : public reduction_2d<min_2d_<T>, T,
                          min_op<typename std::remove_reference_t<T>::Scalar>> {
  using base
      = reduction_2d<min_2d_<T>, T,
                     min_op<typename std::remove_reference_t<T>::Scalar>>;
  using base::arguments_;

 public:
  using op = min_op<typename std::remove_reference_t<T>::Scalar>;
  explicit min_2d_(T&& a)
      : reduction_2d<min_2d_<T>, T, op>(std::forward<T>(a), op::init()) {}
  /**
   * Creates a deep copy of this expression.
   * @return copy of \c *this
   */
  inline auto deep_copy() const {
    auto&& arg_copy = this->template get_arg<0>().deep_copy();
    return min_2d_<std::remove_reference_t<decltype(arg_copy)>>(
        std::move(arg_copy));
  }
};

/**
 * Two dimensional min - reduction of a kernel generator expression.  So as to
 * be efficient two dimensional reductions are only done partially. That means
 * instead of 1 row kernel output will have a few rows that need to be reduced
 * to obtain final result (actually it is 1 result per work group run - roughly
 * 16 times the number of compute units on the OpenCL device). This can be done
 * in a separate kernel or after copying to CPU. Also two dimensional reductions
 * can not be used as arguments to other operations - they can only be
 * evaluated.
 * @tparam T type of input expression
 * @param a the expression to reduce
 * @return min
 */
template <typename T, require_all_kernel_expressions_t<T>* = nullptr>
inline auto min_2d(T&& a) {
  return min_2d_<as_operation_cl_t<T>>(as_operation_cl(std::forward<T>(a)));
}

namespace internal {
template <typename T>
struct is_reduction_2d_impl
    : public std::is_base_of<internal::reduction_2d_base, std::decay_t<T>> {};
template <typename T>
struct is_reduction_2d_impl<calc_if_<true, T>>
    : public std::is_base_of<internal::reduction_2d_base, std::decay_t<T>> {};
}  // namespace internal

/**
 * Check whether a kernel generator expression is a colwise reduction.
 */
template <typename T>
using is_reduction_2d = internal::is_reduction_2d_impl<std::decay_t<T>>;

/** @}*/
}  // namespace math
}  // namespace stan
#endif
#endif
