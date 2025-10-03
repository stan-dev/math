#ifndef STAN_MATH_FWD_FUN_PHI_APPROX_HPP
#define STAN_MATH_FWD_FUN_PHI_APPROX_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/inv_logit.hpp>
#include <stan/math/fwd/fun/pow.hpp>
#include <stan/math/prim/fun/Phi_approx.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return an approximation of the unit normal cumulative
 * distribution function (CDF).
 *
 * @tparam T scalar type of forward-mode autodiff variable
 * argument.
 * @param x argument
 * @return approximate probability random sample is less than or
 * equal to argument
 */
template <typename T>
inline fvar<T> Phi_approx(const fvar<T>& x) {
  return inv_logit(0.07056 * pow(x, 3.0) + 1.5976 * x);
}

}  // namespace math
}  // namespace stan
#endif
