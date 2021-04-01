#ifndef STAN_MATH_PRIM_PROB_POISSON_LOG_RNG_HPP
#define STAN_MATH_PRIM_PROB_POISSON_LOG_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a Poisson random variate with specified log rate parameter
 * using the given random number generator.
 *
 * lambda can be a scalar or a one-dimensional container.
 *
 * @tparam T_rate type of log rate parameter
 * @tparam RNG type of random number generator
 * @param alpha (Sequence of) log rate parameter(s)
 * @param rng random number generator
 * @return (Sequence of) Poisson random variate(s)
 * @throw std::domain_error if alpha is nonfinite
 */
template <typename T_rate, class RNG>
inline typename VectorBuilder<true, int, T_rate>::type poisson_log_rng(
    const T_rate& alpha, RNG& rng) {
  using boost::variate_generator;
  using boost::random::poisson_distribution;
  static const char* function = "poisson_log_rng";
  static const double POISSON_MAX_LOG_RATE = 30 * LOG_TWO;
  const auto& alpha_ref = to_ref(alpha);
  check_finite(function, "Log rate parameter", alpha_ref);
  check_less(function, "Log rate parameter", alpha_ref, POISSON_MAX_LOG_RATE);

  scalar_seq_view<T_rate> alpha_vec(alpha_ref);
  size_t N = stan::math::size(alpha);
  VectorBuilder<true, int, T_rate> output(N);

  for (size_t n = 0; n < N; ++n) {
    variate_generator<RNG&, poisson_distribution<> > poisson_rng(
        rng, poisson_distribution<>(std::exp(alpha_vec[n])));
    output[n] = poisson_rng();
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
