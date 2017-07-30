#ifndef STAN_MATH_PRIM_SCAL_PROB_EXPONENTIAL_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_EXPONENTIAL_RNG_HPP

#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/scalar_seq_view.hpp>
#include <stan/math/prim/scal/meta/max_size.hpp>

namespace stan {
  namespace math {
    /**
     * Return a pseudorandom exponential variate for the given inverse scale
     * using the specified random number generator.
     *
     * beta can be either a scalar or vector.
     *
     * @tparam T_inv_scale Type of beta, can either be a scalar or a vector
     * @tparam RNG class of random number generator
     * @param beta (Sequence of) inverse scale parameter(s)
     * @param rng random number generator
     * @return exponential random variate
     * @throw std::domain_error if beta is nonpositive
     */
    template <typename T_inv_scale, class RNG>
    inline typename VectorBuilder<true, double, T_inv_scale>::type
    exponential_rng(const T_inv_scale& beta,
                    RNG& rng) {
      using boost::variate_generator;
      using boost::exponential_distribution;

      static const char* function("exponential_rng");

      scalar_seq_view<T_inv_scale> beta_vec(beta);
      check_positive_finite(function, "Inverse scale parameter", beta);

      size_t N = length(beta);
      VectorBuilder<true, double, T_inv_scale> output(N);
      for (size_t n = 0; n < N; n++) {
        variate_generator<RNG&, exponential_distribution<> >
          exp_rng(rng, exponential_distribution<>(beta_vec[n]));
        output[n] = exp_rng();
      }

      return output.data();
    }
  }
}
#endif
