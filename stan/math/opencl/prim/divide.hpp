#ifndef STAN_MATH_OPENCL_PRIM_DIVIDE_HPP
#define STAN_MATH_OPENCL_PRIM_DIVIDE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {
/** \ingroup opencl
 * Returns the elementwise division of the kernel generator expression
 *
 * @tparam T_a type of input kernel generator expression a
 * @param a expression to divide
 * @param d scalar to divide by
 * @return the elements of expression a divided by d
 */
template <typename T_a,
          typename = require_all_kernel_expressions_and_none_scalar_t<T_a>>
inline auto divide(T_a&& a, double d) {  // NOLINT
  return elt_divide(std::forward<T_a>(a), d);
}
}  // namespace math
}  // namespace stan
#endif
#endif
