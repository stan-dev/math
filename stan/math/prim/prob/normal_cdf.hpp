#ifndef STAN_MATH_PRIM_PROB_NORMAL_CDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_CDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/prob/normal_lcdf.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Calculates the normal cumulative distribution function for the given
 * variate, location, and scale.
 *
 * \f$\Phi(x) = \frac{1}{\sqrt{2 \pi}} \int_{-\inf}^x e^{-t^2/2} dt\f$.
 *
 * Evaluated as the exponential of the log cdf so the tails and gradients
 * share the log kernel instead of hard cutoffs.
 *
 * @tparam T_y type of y
 * @tparam T_loc type of mean parameter
 * @tparam T_scale type of standard deviation parameter
 * @param y A scalar variate.
 * @param mu The location of the normal distribution.
 * @param sigma The scale of the normal distribution
 * @return The unit normal cdf evaluated at the specified arguments.
 */
template <typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline return_type_t<T_y, T_loc, T_scale> normal_cdf(T_y&& y, T_loc&& mu,
                                                     T_scale&& sigma) {
  return exp(internal::normal_lcdf_impl<false>(
      "normal_cdf", std::forward<T_y>(y), std::forward<T_loc>(mu),
      std::forward<T_scale>(sigma)));
}

}  // namespace math
}  // namespace stan
#endif
