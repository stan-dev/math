#ifndef STAN_MATH_PRIM_PROB_BERNOULLI_LOGIT_RNG_HPP
#define STAN_MATH_PRIM_PROB_BERNOULLI_LOGIT_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a Bernoulli random variate with logit-parameterized chance of success
 * using the specified random number generator.
 *
 * t can be a scalar or a one-dimensional container.
 *
 * @tparam T_t type of logit-parameterized chance of success parameter
 * @tparam RNG type of random number generator
 * @param t (Sequence of) logit-parameterized chance of success parameter(s)
 * @param rng random number generator
 * @return (Sequence of) Bernoulli random variate(s)
 * @throw std::domain_error if logit-parameterized chance of success parameter
 * is not finite
 */
template <typename T_t, class RNG>
inline typename VectorBuilder<true, int, T_t>::type bernoulli_logit_rng(
    const T_t& t, RNG& rng) {
  using boost::bernoulli_distribution;
  using boost::variate_generator;

  check_finite("bernoulli_logit_rng", "Logit transformed probability parameter",
               t);

  scalar_seq_view<T_t> t_vec(t);
  size_t N = size(t);
  VectorBuilder<true, int, T_t> output(N);

  for (size_t n = 0; n < N; ++n) {
    variate_generator<RNG&, bernoulli_distribution<> > bernoulli_rng(
        rng, bernoulli_distribution<>(inv_logit(t_vec[n])));
    output[n] = bernoulli_rng();
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
