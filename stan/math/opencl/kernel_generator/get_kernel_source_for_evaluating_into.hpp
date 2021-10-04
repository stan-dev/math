#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_GET_KERNEL_SOURCE_FOR_EVALUATING_INTO_HPP  // NOLINT
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_GET_KERNEL_SOURCE_FOR_EVALUATING_INTO_HPP  // NOLINT
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator/operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/as_operation_cl.hpp>
#include <stan/math/opencl/kernel_generator/multi_result_kernel.hpp>
#include <CL/opencl.hpp>
#include <string>
#include <set>

namespace stan {
namespace math {

/** \addtogroup opencl_kernel_generator
 *  @{
 */
template <typename Derived, typename Scalar, typename... Args>
template <typename T_lhs>
std::string
operation_cl<Derived, Scalar, Args...>::get_kernel_source_for_evaluating_into(
    const T_lhs& lhs) const {
  static_assert(
      is_kernel_expression<T_lhs>::value,
      "operation_cl::get_kernel_source_for_evaluating_into: left hand "
      "side is not a valid expression!");
  return results(lhs).get_kernel_source_for_evaluating(expressions(derived()));
}
/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
