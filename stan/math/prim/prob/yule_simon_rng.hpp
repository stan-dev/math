#ifndef STAN_MATH_PRIM_PROB_YULE_SIMON_RNG_HPP
#define STAN_MATH_PRIM_PROB_YULE_SIMON_RNG_HPP

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
inline auto yule_simon_rng(const T_alpha &alpha, RNG &rng) {
  using T_alpha_ref = ref_type_t<T_alpha>;
  static constexpr const char *function = "yule_simon_rng";

  T_alpha_ref alpha_ref = alpha;
  check_positive_finite(function, "Shape parameter", alpha_ref);

  using T_w = decltype(exponential_rng(alpha_ref, rng));
  T_w w = exponential_rng(alpha_ref, rng);

  scalar_seq_view<T_w> w_vec(w);
  size_t size_w = stan::math::size(w);

  VectorBuilder<true, int, T_alpha> output(size_w);
  for (size_t n = 0; n < size_w; ++n) {
    double p = stan::math::exp(-w_vec.val(n));
    double odds_ratio_p
        = stan::math::exp(stan::math::log(p) - stan::math::log1m(p));
    output[n] = neg_binomial_rng(1.0, odds_ratio_p, rng) + 1;
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
