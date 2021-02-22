#ifndef STAN_MATH_OPENCL_REV_REVERSE_HPP
#define STAN_MATH_OPENCL_REV_REVERSE_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/block.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Return reversed view into the specified vector or row vector.
 *
 * @tparam T type of expression
 * @param m vector or row vector to reverse
 * @return reversed vector or row vector
 */
template <typename T,
          require_all_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
          require_any_var_t<T>* = nullptr>
inline auto reverse(const T& m) {
  return m.reverse();
}

}  // namespace math
}  // namespace stan

#endif
#endif
