#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_CALC_IF_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_CALC_IF_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl_view.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/opencl/kernel_generator/type_str.hpp>
#include <stan/math/opencl/kernel_generator/name_generator.hpp>
#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/is_kernel_expression.hpp>
#include <string>
#include <type_traits>
#include <set>
#include <utility>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */

/**
 * Represents a calc_if in kernel generator expressions.
 * @tparam T type of the argument
 * @tparam Do_Calculate whether to calculate this expression in a kernel
 */
template <bool Do_Calculate, typename T>
class calc_if_
    : public operation_cl<calc_if_<Do_Calculate, T>,
                          typename std::remove_reference_t<T>::Scalar, T> {
 public:
  using Scalar = typename std::remove_reference_t<T>::Scalar;
  using base = operation_cl<calc_if_<Do_Calculate, T>, Scalar, T>;
  using base::var_name_;

  /**
   * Constructor
   * @param a expression to calc_if
   */
  explicit calc_if_(T&& a) : base(std::forward<T>(a)) {}

  inline kernel_parts generate(const std::string& row_index_name,
                               const std::string& col_index_name,
                               const bool view_handled,
                               const std::string& var_name_arg) const {
    if (Do_Calculate) {
      var_name_ = var_name_arg;
    }
    return {};
  }

  /**
   * Generates kernel code for assigning this expression into result expression.
   * @param[in,out] generated set of (pointer to) already generated operations
   * @param ng name generator for this kernel
   * @param row_index_name row index variable name
   * @param col_index_name column index variable name
   * @param result expression into which result is to be assigned
   * @return part of kernel with code for this and nested expressions
   * @throws std::invalid_argument dimensions of expression and result can not
   * be resized.
   */
  template <typename T_result>
  kernel_parts get_whole_kernel_parts(
      std::set<const operation_cl_base*>& generated, name_generator& ng,
      const std::string& row_index_name, const std::string& col_index_name,
      const T_result& result) const {
    if (Do_Calculate) {
      return this->template get_arg<0>().get_whole_kernel_parts(
          generated, ng, row_index_name, col_index_name, result);
    } else {
      return {};
    }
  }

  /**
   * Sets kernel arguments for nested expressions.
   * @param[in,out] generated set of expressions that already set their kernel
   * arguments
   * @param kernel kernel to set arguments on
   * @param[in,out] arg_num consecutive number of the first argument to set.
   * This is incremented for each argument set by this function.
   */
  inline void set_args(std::set<const operation_cl_base*>& generated,
                       cl::Kernel& kernel, int& arg_num) const {
    if (Do_Calculate) {
      this->template get_arg<0>().set_args(generated, kernel, arg_num);
    }
  }
};

template <bool Do_Calculate, typename T,
          typename = require_all_kernel_expressions_and_none_scalar_t<T>>
inline calc_if_<Do_Calculate, as_operation_cl_t<T>> calc_if(T&& a) {
  return calc_if_<Do_Calculate, as_operation_cl_t<T>>(
      as_operation_cl(std::forward<T>(a)));
}

namespace internal {
template <typename T>
struct is_without_output_impl : std::false_type {};

template <typename T>
struct is_without_output_impl<calc_if_<false, T>> : std::true_type {};
}  // namespace internal

template <typename T>
using is_without_output = internal::is_without_output_impl<std::decay_t<T>>;
/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif  // STAN_MATH_OPENCL_KERNEL_GENERATOR_calc_if_HPP
