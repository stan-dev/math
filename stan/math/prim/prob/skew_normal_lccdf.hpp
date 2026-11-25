#ifndef STAN_MATH_PRIM_PROB_SKEW_NORMAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_SKEW_NORMAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/minus.hpp>
#include <stan/math/prim/prob/skew_normal_lcdf.hpp>

namespace stan {
namespace math {
namespace internal {
constexpr char skew_normal_lccdf_func[] = "skew_normal_lccdf";
}  // namespace internal

/** \ingroup prob_dists
 * Returns the log of the skew normal complementary cumulative distribution
 * function. If X ~ SkewNormal(mu, sigma, alpha), then
 * -X ~ SkewNormal(-mu, sigma, -alpha), so the ccdf at y is the cdf of
 * SkewNormal(-mu, sigma, -alpha) at -y. The upper tail is then computed
 * without the subtraction 1 - cdf.
 *
 * @tparam T_y type of the random variable
 * @tparam T_loc type of the location parameter
 * @tparam T_scale type of the scale parameter
 * @tparam T_shape type of the shape parameter
 * @param y random variable
 * @param mu location parameter
 * @param sigma scale parameter
 * @param alpha shape parameter
 * @return The log of the skew normal ccdf evaluated at the specified
 *   arguments. If given containers, the log of the product of the ccdfs.
 */
template <typename T_y, typename T_loc, typename T_scale, typename T_shape>
inline return_type_t<T_y, T_loc, T_scale, T_shape> skew_normal_lccdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma, const T_shape& alpha) {
  return internal::skew_normal_lcdf_impl<internal::skew_normal_lccdf_func>(
      minus(y), minus(mu), sigma, minus(alpha));
}

}  // namespace math
}  // namespace stan
#endif
