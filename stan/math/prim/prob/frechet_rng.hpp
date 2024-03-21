#ifndef STAN_MATH_PRIM_PROB_FRECHET_RNG_HPP
#define STAN_MATH_PRIM_PROB_FRECHET_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <boost/random/weibull_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {
/** \ingroup prob_dists
 * Return a pseudorandom Frechet variate for the given shape
 * and scale parameters using the specified random number generator.
 *
 * alpha and sigma can each be a scalar or a one-dimensional container. Any
 * non-scalar inputs must be the same size.
 *
 * @tparam T_shape type of shape parameter
 * @tparam T_scale type of scale parameter
 * @tparam RNG type of random number generator
 * @param alpha (Sequence of) positive shape parameter(s)
 * @param sigma (Sequence of) positive scale parameter(s)
 * @param rng random number generator
 * @return (Sequence of) Frechet random variate(s)
 * @throw std::domain_error if alpha is nonpositive or sigma is nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_shape, typename T_scale, class RNG>
inline typename VectorBuilder<true, double, T_shape, T_scale>::type frechet_rng(
    const T_shape& alpha, const T_scale& sigma, RNG& rng) {
  using boost::variate_generator;
  using boost::random::weibull_distribution;
  using T_alpha_ref = ref_type_t<T_shape>;
  using T_sigma_ref = ref_type_t<T_scale>;
  static constexpr const char* function = "frechet_rng";
  check_consistent_sizes(function, "Shape parameter", alpha, "Scale Parameter",
                         sigma);
  T_alpha_ref alpha_ref = alpha;
  T_sigma_ref sigma_ref = sigma;
  check_positive_finite(function, "Shape parameter", alpha_ref);
  check_positive_finite(function, "Scale parameter", sigma_ref);

  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_sigma_ref> sigma_vec(sigma_ref);
  size_t N = max_size(alpha, sigma);
  VectorBuilder<true, double, T_shape, T_scale> output(N);

  for (size_t n = 0; n < N; ++n) {
    variate_generator<RNG&, weibull_distribution<> > weibull_rng(
        rng, weibull_distribution<>(alpha_vec[n], 1.0 / sigma_vec[n]));
    output[n] = 1 / weibull_rng();
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
