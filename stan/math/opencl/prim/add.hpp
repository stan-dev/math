#ifndef STAN_MATH_OPENCL_PRIM_ADD_HPP
#define STAN_MATH_OPENCL_PRIM_ADD_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {
/** \ingroup opencl
 * Computes the sum on kernel generator expressions.
 *
 * @tparam T_a type of input kernel generator expression a
 * @tparam T_b type of input kernel generator expression b
 * @param a first expression
 * @param b second expression
 * @return the sum of the first and second expression
 *
 * @throw <code>std::invalid_argument</code> if the
 *   dimensions dont match if none of tha arguments are a scalar.
 */
template <typename T_a, typename T_b,
          typename = require_all_kernel_expressions_t<T_a, T_b>,
          typename = require_any_not_arithmetic_t<T_a, T_b>>
inline auto add(T_a&& a, T_b&& b) {  // NOLINT
  return std::forward<T_a>(a) + std::forward<T_b>(b);
}
}  // namespace math
}  // namespace stan
#endif
#endif
