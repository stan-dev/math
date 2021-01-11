#ifndef STAN_MATH_PRIM_PROB_VON_MISES_RNG_HPP
#define STAN_MATH_PRIM_PROB_VON_MISES_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>

namespace stan {
namespace math {

/** \ingroup prob_dists
 * Return a von Mises random variate for the given location and concentration
 * using the specified random number generator.
 *
 * mu and kappa can each be a scalar or a vector. Any non-scalar
 * inputs must be the same length.
 *
 * The algorithm used in von_mises_rng is a modified version of the
 * algorithm in:
 *
 * Efficient Simulation of the von Mises Distribution
 * D. J. Best and N. I. Fisher
 * Journal of the Royal Statistical Society. Series C (Applied Statistics),
 * Vol. 28, No. 2 (1979), pp. 152-157
 *
 * For kappa < 1.4e-8, this reduces to a uniform distribution.
 *
 * @tparam T_loc type of location parameter
 * @tparam T_conc type of scale (concentration) parameter
 * @tparam RNG type of random number generator
 *
 * @param mu (Sequence of) location parameter(s)
 * @param kappa (Sequence of) non-negative scale (concentration) parameter(s)
 * @param rng random number generator
 * @return (Sequence of) von Mises random variate(s)
 * @throw std::domain_error if mu is infinite or kappa is negative
 * @throw std::invalid_argument if non-scalar arguments are of different
 * sizes
 */
template <typename T_loc, typename T_conc, class RNG>
inline typename VectorBuilder<true, double, T_loc, T_conc>::type von_mises_rng(
    const T_loc& mu, const T_conc& kappa, RNG& rng) {
  using boost::variate_generator;
  using boost::random::uniform_real_distribution;
  using T_mu_ref = ref_type_t<T_loc>;
  using T_kappa_ref = ref_type_t<T_conc>;
  static const char* function = "von_mises_rng";
  check_consistent_sizes(function, "Location parameter", mu, "Scale parameter",
                         kappa);
  T_mu_ref mu_ref = mu;
  T_kappa_ref kappa_ref = kappa;

  check_finite(function, "Location parameter", mu_ref);
  check_nonnegative(function, "Scale parameter", kappa_ref);
  check_finite(function, "Scale parameter", kappa_ref);

  scalar_seq_view<T_mu_ref> mu_vec(mu_ref);
  scalar_seq_view<T_kappa_ref> kappa_vec(kappa_ref);
  size_t N = max_size(mu, kappa_ref);
  VectorBuilder<true, double, T_loc, T_conc> output(N);

  variate_generator<RNG&, uniform_real_distribution<> > uniform_rng(
      rng, uniform_real_distribution<>(0.0, 1.0));

  for (size_t n = 0; n < N; ++n) {
    // for kappa sufficiently close to zero, it reduces to a
    // circular uniform distribution centered at mu
    if (kappa_vec[n] < 1.4e-8) {
      output[n] = (uniform_rng() - 0.5) * TWO_PI
                  + std::fmod(std::fmod(mu_vec[n], TWO_PI) + TWO_PI, TWO_PI);
      continue;
    }

    double r = 1 + std::pow((1 + 4 * kappa_vec[n] * kappa_vec[n]), 0.5);
    double rho = 0.5 * (r - std::pow(2 * r, 0.5)) / kappa_vec[n];
    double s = 0.5 * (1 + rho * rho) / rho;

    bool done = false;
    double W;
    while (!done) {
      double Z = std::cos(pi() * uniform_rng());
      W = (1 + s * Z) / (s + Z);
      double Y = kappa_vec[n] * (s - W);
      double U2 = uniform_rng();
      done = Y * (2 - Y) - U2 > 0;

      if (!done) {
        done = std::log(Y / U2) + 1 - Y >= 0;
      }
    }

    double U3 = uniform_rng() - 0.5;
    double sign = ((U3 >= 0) - (U3 <= 0));

    //  it's really an fmod() with a positivity constraint
    output[n] = sign * std::acos(W)
                + std::fmod(std::fmod(mu_vec[n], TWO_PI) + TWO_PI, TWO_PI);
  }

  return output.data();
}

}  // namespace math
}  // namespace stan
#endif
