#ifndef STAN_MATH_PRIM_PROB_YULE_SIMON_RNG_HPP
#define STAN_MATH_PRIM_PROB_YULE_SIMON_RNG_HPP

#include <utility>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <stan/math/prim/prob/exponential_rng.hpp>
#include <stan/math/prim/prob/neg_binomial_rng.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a yule-simon random variate with the given shape parameter,
 * using the given random number generator.
 *
 * alpha can be a scalar or a one-dimensional container.
 *
 * @tparam T_alpha type of shape parameter
 * @tparam RNG type of random number generator
 *
 * @param alpha (Sequence of) shape parameter(s)
 * @param rng random number generator
 * @return (Sequence of) yule-simon random variate(s)
 * @throw std::domain_error if alpha is nonpositive
 */
template <typename T_alpha, typename RNG>
inline auto yule_simon_rng(T_alpha&& alpha, RNG& rng) {
  static constexpr const char* function = "yule_simon_rng";
  decltype(auto) alpha_ref = to_ref(std::forward<T_alpha>(alpha));
  check_positive_finite(function, "Shape parameter", alpha_ref);

  auto w = exponential_rng(std::forward<decltype(alpha_ref)>(alpha_ref), rng);
  auto w_arr = as_array_or_scalar(w);
  const auto p = stan::math::exp(-w_arr);
  const auto odds_ratio_p
      = stan::math::exp(stan::math::log(p) - stan::math::log1m(p));

  if constexpr (is_stan_scalar_v<T_alpha>) {
    return neg_binomial_rng(1.0, odds_ratio_p, rng) + 1;
  } else {
    return to_array_1d(
        as_array_or_scalar(neg_binomial_rng(1.0, std::move(odds_ratio_p), rng))
        + 1);
  }
}

}  // namespace math
}  // namespace stan
#endif
