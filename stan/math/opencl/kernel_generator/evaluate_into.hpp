#ifndef STAN_MATH_OPENCL_KERNEL_GENERATOR_EVALUATE_INTO_HPP
#define STAN_MATH_OPENCL_KERNEL_GENERATOR_EVALUATE_INTO_HPP
#ifdef STAN_OPENCL

#include <stan/math/prim/meta/is_kernel_expression.hpp>
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
void operation_cl<Derived, Scalar, Args...>::evaluate_into(T_lhs& lhs) const {
  static_assert(
      is_kernel_expression<T_lhs>::value,
      "operation_cl::evaluate_into: left hand side is not a valid expression!");
  results(lhs) = expressions(derived());
}
/** @}*/
}  // namespace math
}  // namespace stan

#endif
#endif
