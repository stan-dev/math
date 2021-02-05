#ifndef STAN_MATH_PRIM_PROB_BERNOULLI_RNG_HPP
#define STAN_MATH_PRIM_PROB_BERNOULLI_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a Bernoulli random variate with specified chance of success
 * parameter using the specified random number generator.
 *
 * theta can be a scalar or a one-dimensional container.
 *
 * @tparam T_theta type of chance of success parameter
 * @tparam RNG type of random number generator
 * @param theta (Sequence of) chance of success parameter(s)
 * @param rng random number generator
 * @return (Sequence of) Bernoulli random variate(s)
 * @throw std::domain_error if chance of success parameter is less than zero or
 * greater than one.
 */
template <typename T_theta, class RNG>
inline typename VectorBuilder<true, int, T_theta>::type bernoulli_rng(
    const T_theta& theta, RNG& rng) {
  using boost::bernoulli_distribution;
  using boost::variate_generator;
  static const char* function = "bernoulli_rng";
  ref_type_t<T_theta> theta_ref = theta;
  check_bounded(function, "Probability parameter", value_of(theta_ref), 0.0,
                1.0);

  scalar_seq_view<T_theta> theta_vec(theta_ref);
  size_t N = stan::math::size(theta);
  VectorBuilder<true, int, T_theta> output(N);

  for (size_t n = 0; n < N; ++n) {
    variate_generator<RNG&, bernoulli_distribution<> > bernoulli_rng(
        rng, bernoulli_distribution<>(theta_vec[n]));
    output[n] = bernoulli_rng();
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
