#ifndef STAN_MATH_PRIM_SCAL_PROB_NORMAL_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NORMAL_RNG_HPP

#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/max_size.hpp>

namespace stan {
  namespace math {
    /**
     * Return a pseudorandom Normal variate for the given location and scale
     * using the specified random number generator.
     *
     * mu and sigma can be mixes of either scalars or vector types. If mu and
     * sigma are vector types, they both need to be the same size.
     *
     * @tparam T_loc Type of mu, can either be scalar or vector
     * @tparam T_scale Type of sigma, can either be scalar or vector
     * @tparam RNG type of random number generator
     * @param mu location parameter
     * @param sigma positive scale parameter
     * @param rng random number generator
     * @return Normal random variate
     * @throw std::domain_error if mu is infinite or sigma is nonpositive
     * @throw std::invalid_argument if mu and sigma are vector types of
     * different sizes
     */
    template <typename T_loc, typename T_scale, class RNG>
    inline typename VectorBuilder<true, double, T_loc, T_scale>::type
    normal_rng(const T_loc &mu,
               const T_scale &sigma,
               RNG& rng) {
      using boost::variate_generator;
      using boost::normal_distribution;

      static const char* function("normal_rng");

      scalar_seq_view<T_loc> mu_vec(mu);
      scalar_seq_view<T_scale> sigma_vec(sigma);

      check_finite(function, "Location parameter", mu);
      check_positive_finite(function, "Scale parameter", sigma);
      check_consistent_sizes(function,
                             "Location parameter", mu,
                             "Scale Parameter", sigma);

      size_t N = max_size(mu, sigma);

      VectorBuilder<true, double, T_loc, T_scale> output(N);

      for (size_t n = 0; n < N; n++) {
        variate_generator<RNG&, normal_distribution<> >
          norm_rng(rng, normal_distribution<>(mu_vec[n], sigma_vec[n]));
        output[n] = norm_rng();
      }

      return output.data();
    }
  }
}
#endif
