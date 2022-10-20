#ifndef STAN_MATH_PRIM_PROB_BETA_RNG_HPP
#define STAN_MATH_PRIM_PROB_BETA_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log_sum_exp.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a Beta random variate with the supplied success and failure
 * parameters using the given random number generator.
 *
 * alpha and beta can each be a scalar or a one-dimensional container. Any
 * non-scalar inputs must be the same size.
 *
 * @tparam T_shape1 type of success parameter
 * @tparam T_shape2 type of failure parameter
 * @tparam RNG type of random number generator
 * @param alpha (Sequence of) positive finite success parameter(s)
 * @param beta (Sequence of) positive finite failure parameter(s)
 * @param rng random number generator
 * @return (Sequence of) beta random variate(s)
 * @throw std::domain_error if alpha or beta are nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_shape1, typename T_shape2, class RNG>
inline typename VectorBuilder<true, double, T_shape1, T_shape2>::type beta_rng(
    const T_shape1 &alpha, const T_shape2 &beta, RNG &rng) {
  using boost::variate_generator;
  using boost::random::gamma_distribution;
  using boost::random::uniform_real_distribution;
  using T_alpha_ref = ref_type_t<T_shape1>;
  using T_beta_ref = ref_type_t<T_shape2>;
  static constexpr const char *function = "beta_rng";
  check_consistent_sizes(function, "First shape parameter", alpha,
                         "Second shape Parameter", beta);
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_positive_finite(function, "First shape parameter", alpha_ref);
  check_positive_finite(function, "Second shape parameter", beta_ref);

  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t N = max_size(alpha, beta);
  VectorBuilder<true, double, T_shape1, T_shape2> output(N);

  variate_generator<RNG &, uniform_real_distribution<>> uniform_rng(
      rng, uniform_real_distribution<>(0.0, 1.0));
  for (size_t n = 0; n < N; ++n) {
    // If alpha and beta are large, trust the usual ratio of gammas
    // method for generating beta random variables. If any parameter
    // is small, work in log space and use Marsaglia and Tsang's trick
    if (alpha_vec[n] > 1.0 && beta_vec[n] > 1.0) {
      variate_generator<RNG &, gamma_distribution<>> rng_gamma_alpha(
          rng, gamma_distribution<>(alpha_vec[n], 1.0));
      variate_generator<RNG &, gamma_distribution<>> rng_gamma_beta(
          rng, gamma_distribution<>(beta_vec[n], 1.0));
      double a = rng_gamma_alpha();
      double b = rng_gamma_beta();
      output[n] = a / (a + b);
    } else {
      variate_generator<RNG &, gamma_distribution<>> rng_gamma_alpha(
          rng, gamma_distribution<>(alpha_vec[n] + 1, 1.0));
      variate_generator<RNG &, gamma_distribution<>> rng_gamma_beta(
          rng, gamma_distribution<>(beta_vec[n] + 1, 1.0));
      double log_a = std::log(uniform_rng()) / alpha_vec[n]
                     + std::log(rng_gamma_alpha());
      double log_b
          = std::log(uniform_rng()) / beta_vec[n] + std::log(rng_gamma_beta());
      double log_sum = log_sum_exp(log_a, log_b);
      output[n] = std::exp(log_a - log_sum);
    }
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
