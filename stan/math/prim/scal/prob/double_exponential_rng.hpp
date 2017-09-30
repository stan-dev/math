#ifndef STAN_MATH_PRIM_SCAL_PROB_DOUBLE_EXPONENTIAL_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_DOUBLE_EXPONENTIAL_RNG_HPP

#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/meta/max_size.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <string>

namespace stan {
  namespace math {
    /**
     * Return a pseudorandom double exponential variate with the given location
     * and scale using the specified random number generator.
     *
     * mu and sigma can each be either scalars or vectors. All vector inputs
     * must be the same length.
     *
     * @tparam T_loc Type of location parameter
     * @tparam T_scale Type of scale parameter
     * @tparam RNG class of random number generator
     * @param mu (Sequence of) location parameter(s)
     * @param sigma (Sequence of) scale parameter(s)
     * @param rng random number generator
     * @return double exponential random variate
     * @throw std::domain_error if mu is infinite or sigma is nonpositive
     * @throw std::invalid_argument if vector arguments are not the same length
     */
    template <typename T_loc, typename T_scale, class RNG>
    inline typename VectorBuilder<true, double, T_loc, T_scale>::type
    double_exponential_rng(const T_loc& mu,
                           const T_scale& sigma,
                           RNG& rng) {
      static const std::string function = "double_exponential_rng";

      using boost::variate_generator;
      using boost::random::uniform_01;
      using std::abs;

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
        variate_generator<RNG&, uniform_01<> > rng_unit_01(rng, uniform_01<>());
        double a = 0;
        double laplaceRN = rng_unit_01();
        if (0.5 - laplaceRN > 0)
          a = 1.0;
        else if (0.5 - laplaceRN < 0)
          a = -1.0;
        output[n] = mu_vec[n] -
          sigma_vec[n] * a * log1m(2 * abs(0.5 - laplaceRN));
      }

      return output.data();
    }
  }
}
#endif
