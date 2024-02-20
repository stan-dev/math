#ifndef STAN_MATH_PRIM_PROB_SCALED_INV_CHI_SQUARE_RNG_HPP
#define STAN_MATH_PRIM_PROB_SCALED_INV_CHI_SQUARE_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <boost/random/chi_squared_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a scaled chi square random variate for the given
 * number of degrees of freedom and scale using the specified random
 * number generator.
 *
 * nu and sigma can each be a scalar or a one-dimensional container. Any
 * non-scalar inputs must be the same size.
 *
 * @tparam T_deg type of degrees of freedom parameter
 * @tparam T_scale type of scale parameter
 * @tparam RNG type of random number generator
 * @param nu (Sequence of) positive degrees of freedom parameter(s)
 * @param s (Sequence of) positive scale parameter(s)
 * @param rng random number generator
 * @return (Sequence of) scaled chi square random variate(s)
 * @throw std::domain_error if nu or sigma are nonpositive
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_deg, typename T_scale, class RNG>
inline typename VectorBuilder<true, double, T_deg, T_scale>::type
scaled_inv_chi_square_rng(const T_deg& nu, const T_scale& s, RNG& rng) {
  using boost::variate_generator;
  using boost::random::chi_squared_distribution;
  static constexpr const char* function = "scaled_inv_chi_square_rng";
  check_consistent_sizes(function, "Location parameter", nu, "Scale Parameter",
                         s);
  const auto& nu_ref = to_ref(nu);
  const auto& s_ref = to_ref(s);
  check_positive_finite(function, "Degrees of freedom parameter", nu_ref);
  check_positive_finite(function, "Scale parameter", s_ref);

  scalar_seq_view<T_deg> nu_vec(nu_ref);
  scalar_seq_view<T_scale> s_vec(s_ref);
  size_t N = max_size(nu, s);
  VectorBuilder<true, double, T_deg, T_scale> output(N);

  for (size_t n = 0; n < N; ++n) {
    variate_generator<RNG&, chi_squared_distribution<> > chi_square_rng(
        rng, chi_squared_distribution<>(nu_vec[n]));
    output[n] = nu_vec[n] * s_vec[n] * s_vec[n] / chi_square_rng();
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
