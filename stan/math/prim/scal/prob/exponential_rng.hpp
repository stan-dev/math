#ifndef STAN_MATH_PRIM_SCAL_PROB_EXPONENTIAL_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_EXPONENTIAL_RNG_HPP

#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {
/**
 * Return a pseudorandom exponential variate with inverse scale beta
 * using the specified random number generator.
 *
 * beta can be a scalar, a std::vector, an Eigen::Vector, or an
 * Eigen::RowVector.
 *
 * @tparam T_inv Type of inverse scale parameter
 * @tparam RNG class of random number generator
 * @param beta (Sequence of) positive inverse scale parameter(s)
 * @param rng random number generator
 * @return exponential random variate
 * @throw std::domain_error if beta is nonpositive
 */
template <typename T_inv, class RNG>
inline typename VectorBuilder<true, double, T_inv>::type exponential_rng(
    const T_inv& beta, RNG& rng) {
  using boost::variate_generator;
  using boost::exponential_distribution;
  static const char* function = "exponential_rng";

  check_positive_finite(function, "Inverse scale parameter", beta);

  scalar_seq_view<T_inv> beta_vec(beta);
  size_t N = length(beta);
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

}  // namespace math
}  // namespace stan
#endif
