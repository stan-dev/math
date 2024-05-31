#ifndef STAN_MATH_OPENCL_PRIM_DOT_SELF_HPP
#define STAN_MATH_OPENCL_PRIM_DOT_SELF_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/prim/sum.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/check_matching_sizes.hpp>

namespace stan {
namespace math {

/**
 * Returns squared norm of a vector or matrix. For vectors that equals the dot
 * product of the specified vector with itself.
 *
 * @tparam T type of the vector
 * @param a Vector.
 */
template <typename T,
          require_all_kernel_expressions_and_none_scalar_t<T>* = nullptr>
inline auto dot_self(const T& a) {
  return sum(elt_multiply(a, a));
}

}  // namespace math
}  // namespace stan

#endif
#endif
