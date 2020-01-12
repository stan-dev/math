#ifndef STAN_MATH_PRIM_PROB_POISSON_RNG_HPP
#define STAN_MATH_PRIM_PROB_POISSON_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a Poisson random variate with specified rate parameter
 * using the given random number generator.
 *
 * lambda can be a scalar or a one-dimensional container.
 *
 * @tparam T_rate type of rate parameter
 * @tparam RNG type of random number generator
 * @param lambda (Sequence of) rate parameter(s)
 * @param rng random number generator
 * @return (Sequence of) Poisson random variate(s)
 * @throw std::domain_error if lambda is nonpositive
 */
template <typename T_rate, class RNG>
inline typename VectorBuilder<true, int, T_rate>::type poisson_rng(
    const T_rate& lambda, RNG& rng) {
  using boost::random::poisson_distribution;
  using boost::variate_generator;

  static const char* function = "poisson_rng";

  check_not_nan(function, "Rate parameter", lambda);
  check_positive(function, "Rate parameter", lambda);
  check_less(function, "Rate parameter", lambda, POISSON_MAX_RATE);

  scalar_seq_view<T_rate> lambda_vec(lambda);
  size_t N = size(lambda);
  VectorBuilder<true, int, T_rate> output(N);

  for (size_t n = 0; n < N; ++n) {
    variate_generator<RNG&, poisson_distribution<> > poisson_rng(
        rng, poisson_distribution<>(lambda_vec[n]));
    output[n] = poisson_rng();
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
