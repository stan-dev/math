#ifndef STAN_MATH_PRIM_SCAL_FUN_LOCSCALE_FREE_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LOCSCALE_FREE_HPP

#include <stan/math/prim/scal/fun/identity_free.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <boost/math/tools/promotion.hpp>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * Return the unconstrained scalar that transforms to the
 * specified lower- and upper-bounded scalar given the specified
 * bounds.
 *
 * <p>The transfrom in <code>lub_constrain(T, double, double)</code>,
 * is reversed by a transformed and scaled logit,
 *
 * <p>\f$f^{-1}(y) = \mbox{logit}(\frac{y - L}{U - L})\f$
 *
 * where \f$U\f$ and \f$L\f$ are the lower and upper bounds.
 *
 * <p>If the lower bound is negative infinity and upper bound finite,
 * this function reduces to <code>ub_free(y, ub)</code>.  If
 * the upper bound is positive infinity and the lower bound
 * finite, this function reduces to
 * <code>lb_free(x, lb)</code>.  If the upper bound is
 * positive infinity and the lower bound negative infinity,
 * this function reduces to <code>identity_free(y)</code>.
 *
 * @tparam T type of scalar
 * @tparam M type of mean
 * @tparam S type of scale
 * @param y constrained value
 * @param[in] mu location of constrained output
 * @param[in] sigma scale of constrained output
 * @return the free scalar that transforms to the input scalar
 *   given the location and scale
 * @throw std::domain_error if sigma <= 0
 */
template <typename T, typename M, typename S>
inline typename boost::math::tools::promote_args<T, M, S>::type locscale_free(
    const T& y, const M& mu, const S& sigma) {
  if (sigma == 1)
    if (mu == 0)
      return identity_free(y);
  check_positive_finite("locscale_constrain", "scale", sigma);
  return (y - mu) / sigma;
}

}  // namespace math
}  // namespace stan
#endif
