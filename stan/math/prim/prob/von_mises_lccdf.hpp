#ifndef STAN_MATH_PRIM_PROB_VON_MISES_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_VON_MISES_LCCDF_HPP

#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/prob.hpp>
#include <stan/math/prim/fun/modified_bessel_first_kind.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Calculates the log of the complement of the cumulative distribution
 * function of the von Mises distribution:
 *
 * \f$VonMisesLCCDF(x, \mu, \kappa) = \log ( 1.0 - \frac{1}{2\pi I_0(\kappa)}
 * \int_{-\pi}^x\f$ \f$e^{\kappa cos(t - \mu)} dt )\f$
 *
 * where
 *
 * \f$x \in [-\pi, \pi]\f$, \f$\mu \in \mathbb{R}\f$,
 * and \f$\kappa \in \mathbb{R}^+\f$.
 *
 * @param x A scalar variate on the interval \f$[-pi, pi]\f$
 * @param mu The mean of the distribution
 * @param k The inverse scale of the distriubtion
 * @return The log of the von Mises cdf evaluated at the specified arguments
 * @tparam T_x Type of x
 * @tparam T_mu Type of mean parameter
 * @tparam T_k Type of inverse scale parameter
 */
template <typename T_x, typename T_mu, typename T_k>
inline return_type_t<T_x, T_mu, T_k> von_mises_lccdf(const T_x& x,
                                                     const T_mu& mu,
                                                     const T_k& k) {
  return log1m(von_mises_cdf(x, mu, k));
}

}  // namespace math
}  // namespace stan

#endif
