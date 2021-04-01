#ifndef STAN_MATH_PRIM_PROB_EXPONENTIAL_RNG_HPP
#define STAN_MATH_PRIM_PROB_EXPONENTIAL_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a exponential random variate with inverse scale beta
 * using the specified random number generator.
 *
 * beta can be a scalar or a one-dimensional container.
 *
 * @tparam T_inv Type of inverse scale parameter
 * @tparam RNG class of random number generator
 * @param beta (Sequence of) positive inverse scale parameter(s)
 * @param rng random number generator
 * @return (Sequence of) exponential random variate(s)
 * @throw std::domain_error if beta is nonpositive
 */
template <typename T_inv, class RNG>
inline typename VectorBuilder<true, double, T_inv>::type exponential_rng(
    const T_inv& beta, RNG& rng) {
  using boost::exponential_distribution;
  using boost::variate_generator;
  static const char* function = "exponential_rng";
  using T_beta_ref = ref_type_t<T_inv>;
  T_beta_ref beta_ref = beta;
  check_positive_finite(function, "Inverse scale parameter", beta_ref);

  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t N = stan::math::size(beta);
  VectorBuilder<true, double, T_inv> output(N);

  for (size_t n = 0; n < N; ++n) {
    variate_generator<RNG&, exponential_distribution<> > exp_rng(
        rng, exponential_distribution<>(beta_vec[n]));
    output[n] = exp_rng();
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
