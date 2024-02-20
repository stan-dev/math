#ifndef STAN_MATH_PRIM_PROB_LOGLOGISTIC_RNG_HPP
#define STAN_MATH_PRIM_PROB_LOGLOGISTIC_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a loglogistic random variate for the given scale and
 * shape parameters using the specified random number generator.
 *
 * alpha and beta can each be a scalar or a one-dimensional container. Any
 * non-scalar inputs must be the same size.
 *
 * @tparam T_scale type of scale parameter
 * @tparam T_shape type of shape parameter
 * @tparam RNG type of random number generator
 * @param alpha (Sequence of) positive scale parameter(s)
 * @param beta (Sequence of) positive shape parameter(s)
 * @param rng random number generator
 * @return (Sequence of) loglogistic random variate(s)
 * @throw std::domain_error if alpha or beta are nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_scale, typename T_shape, class RNG>
inline typename VectorBuilder<true, double, T_scale, T_shape>::type
loglogistic_rng(const T_scale& alpha, const T_shape& beta, RNG& rng) {
  using boost::uniform_01;
  using boost::variate_generator;
  using std::pow;
  using T_alpha_ref = ref_type_t<T_scale>;
  using T_beta_ref = ref_type_t<T_shape>;
  static constexpr const char* function = "loglogistic_rng";
  check_consistent_sizes(function, "Scale parameter", alpha, "Shape Parameter",
                         beta);
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;
  check_positive_finite(function, "Scale parameter", alpha_ref);
  check_positive_finite(function, "Shape parameter", beta_ref);

  scalar_seq_view<T_alpha_ref> alpha_vec(alpha_ref);
  scalar_seq_view<T_beta_ref> beta_vec(beta_ref);
  size_t N = max_size(alpha, beta);
  VectorBuilder<true, double, T_scale, T_shape> output(N);

  for (size_t n = 0; n < N; ++n) {
    variate_generator<RNG&, uniform_01<> > uniform01_rng(rng, uniform_01<>());
    const double tmp = uniform01_rng();
    output[n] = alpha_vec[n] * pow(tmp / (1 - tmp), 1 / beta_vec[n]);
  }
  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
