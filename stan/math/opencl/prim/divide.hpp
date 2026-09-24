#ifndef STAN_MATH_OPENCL_PRIM_DIVIDE_HPP
#define STAN_MATH_OPENCL_PRIM_DIVIDE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/scalar_cl.hpp>

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

/** \ingroup opencl
 * Returns the elements of the first argument divided by the device scalar.
 * @tparam T_a type of the expression
 * @tparam T_d type of the device scalar
 * @param a expression to divide
 * @param d device scalar to divide by
 * @return the elements of expression a divided by d
 */
template <typename T_a, typename T_d,
          require_all_kernel_expressions_and_none_scalar_t<T_a>* = nullptr,
          require_prim_scalar_cl_t<T_d>* = nullptr>
inline auto divide(T_a&& a, T_d&& d) {  // NOLINT
  return elt_divide(std::forward<T_a>(a), std::forward<T_d>(d));
}
}  // namespace math
}  // namespace stan
#endif
#endif
