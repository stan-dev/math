#ifndef STAN_MATH_PRIM_PROB_BINOMIAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_BINOMIAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <boost/random/binomial_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a pseudorandom binomial random variable for the given population
 * size and chance of success parameters using the specified random number
 * generator.
 *
 * beta can be a scalar or a one-dimensional container.
 *
 * @tparam T_N Type of population size parameter
 * @tparam T_theta Type of change of success parameter
 * @tparam RNG class of rng
 * @param N (Sequence of) population size parameter(s)
 * @param theta (Sequence of) chance of success parameter(s)
 * @param rng random number generator
 * @return (Sequence of) binomial random variate(s)
 * @throw std::domain_error if N is negative
 * @throw std::domain_error if theta is not a valid probability
 */
template <typename T_N, typename T_theta, class RNG>
inline typename VectorBuilder<true, int, T_N, T_theta>::type binomial_rng(
    const T_N& N, const T_theta& theta, RNG& rng) {
  using boost::binomial_distribution;
  using boost::variate_generator;
  using T_N_ref = ref_type_t<T_N>;
  using T_theta_ref = ref_type_t<T_theta>;
  static constexpr const char* function = "binomial_rng";
  check_consistent_sizes(function, "Population size parameter", N,
                         "Probability Parameter", theta);
  T_N_ref N_ref = N;
  T_theta_ref theta_ref = theta;
  check_nonnegative(function, "Population size parameter", N_ref);
  check_bounded(function, "Probability parameter", value_of(theta_ref), 0.0,
                1.0);

  scalar_seq_view<T_N_ref> N_vec(N_ref);
  scalar_seq_view<T_theta_ref> theta_vec(theta_ref);
  size_t M = max_size(N, theta);
  VectorBuilder<true, int, T_N, T_theta> output(M);

  for (size_t m = 0; m < M; ++m) {
    variate_generator<RNG&, binomial_distribution<> > binomial_rng(
        rng, binomial_distribution<>(N_vec[m], theta_vec[m]));

    output[m] = binomial_rng();
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
