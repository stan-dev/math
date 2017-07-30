#ifndef STAN_MATH_PRIM_SCAL_PROB_SCALED_INV_CHI_SQUARE_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_SCALED_INV_CHI_SQUARE_RNG_HPP

#include <boost/random/chi_squared_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/max_size.hpp>

namespace stan {
  namespace math {
    /**
     * Return a pseudorandom scaled inverse chi squared variate with the given
     * degrees of freedom and scale using the specified random number generator.
     *
     * nu and sigma can each be either scalars or vectors. Any vector inputs
     * must be the same length.
     *
     * @tparam T_deg Type of degrees of freedom parameter
     * @tparam T_scale Type of scale parameter
     * @tparam RNG type of random number generator
     * @param nu (Sequence of) number of degrees of freedom parameter(s)
     * @param sigma (Sequence of) scale parameter(s)
     * @param rng random number generator
     * @return scaled inverse chi squared random variate
     * @throw std::domain_error if nu or sigma are nonpositive
     * @throw std::invalid_argument if vector arguments are not the same length
     */
    template <typename T_deg, typename T_scale, class RNG>
    inline typename VectorBuilder<true, double, T_deg, T_scale>::type
    scaled_inv_chi_square_rng(const T_deg& nu,
                              const T_scale& sigma,
                              RNG& rng) {
      using boost::variate_generator;
      using boost::random::chi_squared_distribution;

      static const char* function("scaled_inv_chi_square_rng");

      scalar_seq_view<T_deg> nu_vec(nu);
      scalar_seq_view<T_scale> sigma_vec(sigma);

      check_positive_finite(function, "Degrees of freedom parameter", nu);
      check_positive_finite(function, "Scale parameter", sigma);
      check_consistent_sizes(function,
                             "Degrees of freedom parameter", nu,
                             "Scale Parameter", sigma);

      size_t N = max_size(nu, sigma);

      VectorBuilder<true, double, T_deg, T_scale> output(N);

      for (size_t n = 0; n < N; n++) {
        variate_generator<RNG&, chi_squared_distribution<> >
          chi_square_rng(rng, chi_squared_distribution<>(nu_vec[n]));
        output[n] = nu_vec[n] * sigma_vec[n] / chi_square_rng();
      }

      return output.data();
    }

  }
}
#endif
