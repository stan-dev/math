#ifndef STAN_MATH_PRIM_SCAL_PROB_SKEW_NORMAL_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_SKEW_NORMAL_RNG_HPP

#include <boost/random/variate_generator.hpp>
#include <boost/math/distributions.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/prob/uniform_rng.hpp>
#include <string>

namespace stan {
  namespace math {
    /**
     * Return a pseudorandom Skew-normal variate for the given location, scale,
     * and shape using the specified random number generator.
     *
     * mu, alpha, and sigma can each be either scalars or vectors. All vector
     * inputs must be the same length.
     *
     * @tparam T_loc Type of location parameter
     * @tparam T_scale Type of scale parameter
     * @tparam T_shape Type of shape parameter
     * @tparam RNG type of random number generator
     * @param mu (Sequence of) location parameter(s)
     * @param sigma (Sequence of) scale parameter(s)
     * @param alpha (Sequence of) shape parameter(s)
     * @param rng random number generator
     * @return Skew-normal random variate
     * @throw std::domain_error if mu is infinite, alpha is infinite, or sigma
     * is nonpositive
     * @throw std::invalid_argument if vector arguments are not the same length
     */
    template <typename T_loc, typename T_scale, typename T_shape, class RNG>
    inline typename VectorBuilder<true, double, T_loc, T_scale, T_shape>::type
    skew_normal_rng(const T_loc& mu,
                    const T_scale& sigma,
                    const T_shape& alpha,
                    RNG& rng) {
      static const std::string function = "skew_normal_rng";

      scalar_seq_view<T_loc> mu_vec(mu);
      scalar_seq_view<T_scale> sigma_vec(sigma);
      scalar_seq_view<T_shape> alpha_vec(alpha);

      check_finite(function, "Location parameter", mu);
      check_finite(function, "Shape parameter", alpha);
      check_positive_finite(function, "Scale parameter", sigma);
      check_consistent_sizes(function,
                             "Location parameter", mu,
                             "Scale Parameter", sigma,
                             "Shape Parameter", alpha);

      size_t N = max_size(mu, sigma, alpha);

      VectorBuilder<true, double, T_loc, T_scale, T_shape> output(N);

      for (size_t n = 0; n < N; n++) {
        boost::math::skew_normal_distribution<> dist(mu_vec[n], sigma_vec[n],
                                                     alpha_vec[n]);
        output[n] = quantile(dist, uniform_rng(0.0, 1.0, rng));
      }

      return output.data();
    }
  }
}
#endif
