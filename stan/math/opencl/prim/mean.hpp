#ifndef STAN_MATH_OPENCL_PRIM_MEAN_HPP
#define STAN_MATH_OPENCL_PRIM_MEAN_HPP
#ifdef STAN_OPENCL

#include <stan/math/opencl/prim/sum.hpp>

namespace stan {
namespace math {

/**
 * Returns the unary plus of the input.
 *
 * @param M input kernel expression
 * @return result of unary plus of the input.
 */
template <typename T,
          require_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr>
scalar_type_t<T> mean(const T& m) {
  check_nonzero_size("mean", "m", m);
  return sum(m)/m.size();
}

}  // namespace math
}  // namespace stan

#endif
#endif
