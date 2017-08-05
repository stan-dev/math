#ifndef STAN_MATH_PRIM_SCAL_PROB_WEIBULL_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_WEIBULL_RNG_HPP

#include <boost/random/weibull_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/max_size.hpp>
#include <string>

namespace stan {
  namespace math {
    /**
     * Return a pseudorandom Weibull variate with given shape and scale using
     * the specified random number generator.
     *
     * alpha and sigma can each be either scalars or vectors. Any vector inputs
     * must be the same length.
     *
     * @tparam T_shape Type of shape parameter
     * @tparam T_scale Type of scale parameter
     * @tparam RNG class of random number generator
     * @param alpha (Sequence of) shape parameter(s)
     * @param sigma (Sequence of) scale parameter(s)
     * @param rng random number generator
     * @return Weibull random variate
     * @throw std::domain_error if alpha or sigma are nonpositive
     * @throw std::invalid_argument if vector arguments are not the same length
     */
    template <typename T_shape, typename T_scale, class RNG>
    inline typename VectorBuilder<true, double, T_shape, T_scale>::type
    weibull_rng(const T_shape& alpha,
                const T_scale& sigma,
                RNG& rng) {
      using boost::variate_generator;
      using boost::random::weibull_distribution;

      static const std::string function = "weibull_rng";

      scalar_seq_view<T_shape> alpha_vec(alpha);
      scalar_seq_view<T_scale> sigma_vec(sigma);

      check_positive_finite(function, "Shape parameter", alpha);
      check_positive_finite(function, "Scale parameter", sigma);
      check_consistent_sizes(function,
                             "Shape parameter", alpha,
                             "Scale parameter", sigma);

      size_t N = max_size(alpha, sigma);

      VectorBuilder<true, double, T_shape, T_scale> output(N);

      for (size_t n = 0; n < N; n++) {
        variate_generator<RNG&, weibull_distribution<> >
          weibull_rng(rng, weibull_distribution<>(alpha_vec[n], sigma_vec[n]));
        output[n] = weibull_rng();
      }

      return output.data();
    }
  }
}
#endif
