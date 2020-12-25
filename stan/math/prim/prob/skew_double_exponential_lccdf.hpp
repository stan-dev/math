#ifndef STAN_MATH_PRIM_PROB_SKEW_DOUBLE_EXPONENTIAL_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_SKEW_DOUBLE_EXPONENTIAL_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Returns the skew double exponential log complementary cumulative density
 * function. Given containers of matching sizes, returns the log sum of
 * probabilities.
 *
 * @tparam T_y type of real parameter.
 * @tparam T_loc type of location parameter.
 * @tparam T_scale type of scale parameter.
 * @tparam T_skewness type of skewness parameter.
 * @param y real parameter
 * @param mu location parameter
 * @param sigma scale parameter
 * @param tau skewness parameter
 * @return log probability or log sum of probabilities
 * @throw std::domain_error if mu is infinite or sigma is nonpositive or tau is
 *  not bound between 0.0 and 1.0
 * @throw std::invalid_argument if container sizes mismatch
 */
template <typename T_y, typename T_loc, typename T_scale, typename T_skewness>
return_type_t<T_y, T_loc, T_scale, T_skewness> skew_double_exponential_lccdf(
    const T_y& y, const T_loc& mu, const T_scale& sigma,
    const T_skewness& tau) {
  return log1m_exp(skew_double_exponential_lcdf(y, mu, sigma, tau));
}
}  // namespace math
}  // namespace stan
#endif
