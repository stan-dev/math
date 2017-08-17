#ifndef STAN_MATH_PRIM_SCAL_PROB_GUMBEL_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_GUMBEL_RNG_HPP

#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/meta/max_size.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>
#include <string>

namespace stan {
  namespace math {
    /**
     * Return a pseudorandom Gumbel variate with the given location and scale
     * using the specified random number generator.
     *
     * mu and sigma can each be either scalars or vectors. All vector inputs
     * must be the same length.
     *
     * @tparam T_loc Type of location parameter
     * @tparam T_scale Type of scale parameter
     * @tparam RNG type of random number generator
     * @param mu (Sequence of) location parameter(s)
     * @param beta (Sequence of) scale parameter(s)
     * @param rng random number generator
     * @return Gumbel random variate
     * @throw std::domain_error if mu is infinite or beta is nonpositive.
     * @throw std::invalid_argument if vector arguments are not the same length
     */
    template <typename T_loc, typename T_scale, class RNG>
    inline typename VectorBuilder<true, double, T_loc, T_scale>::type
    gumbel_rng(const T_loc &mu,
               const T_scale &beta,
               RNG& rng) {
      using boost::variate_generator;
      using boost::uniform_01;

      static const std::string function = "gumbel_rng";

      check_finite(function, "Location parameter", mu);
      check_positive(function, "Scale parameter", beta);
      check_consistent_sizes(function,
                             "Location parameter", mu,
                             "Scale Parameter", beta);

      scalar_seq_view<T_loc> mu_vec(mu);
      scalar_seq_view<T_scale> beta_vec(beta);

      size_t N = max_size(mu, beta);

      VectorBuilder<true, double, T_loc, T_scale> output(N);

      for (size_t n = 0; n < N; n++) {
        variate_generator<RNG&, uniform_01<> >
          uniform01_rng(rng, uniform_01<>());
        output[n] = mu_vec[n] -
          beta_vec[n] * std::log(-std::log(uniform01_rng()));
      }

      return output.data();
    }
  }
}
#endif
