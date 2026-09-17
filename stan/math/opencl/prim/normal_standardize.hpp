#ifndef STAN_MATH_OPENCL_PRIM_NORMAL_STANDARDIZE_HPP
#define STAN_MATH_OPENCL_PRIM_NORMAL_STANDARDIZE_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/kernel_generator.hpp>

namespace stan {
namespace math {
namespace internal {

/** Standardize finite inputs even when y-mu overflows. */
template <typename T_y, typename T_mu, typename T_sigma>
inline auto normal_standardize_cl(const T_y& y, const T_mu& mu,
                                  const T_sigma& sigma) {
  return select(isinf(y - mu) && isfinite(y) && isfinite(mu),
                elt_divide(y, sigma) - elt_divide(mu, sigma),
                elt_divide(y - mu, sigma));
}

}  // namespace internal
}  // namespace math
}  // namespace stan
#endif
#endif
