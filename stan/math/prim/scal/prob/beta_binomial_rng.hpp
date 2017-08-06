#ifndef STAN_MATH_PRIM_SCAL_PROB_BETA_BINOMIAL_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BETA_BINOMIAL_RNG_HPP

#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/prob/binomial_rng.hpp>
#include <stan/math/prim/scal/prob/beta_rng.hpp>
#include <string>

namespace stan {
  namespace math {

    /**
     * Return a pseudorandom Beta-Binomial draw with specified population size,
     * prior success, and prior failure parameter using the specified random
     * number generator.
     *
     * @tparam RNG type of random number generator
     * @param N population size parameter
     * @param alpha success parameter
     * @param beta failure parameter
     * @param rng random number generator
     * @return Beta-Binomial random variate
     * @throw std::domain_error if N, alpha, or beta is negative.
     */
    template <class RNG>
    inline int
    beta_binomial_rng(int N,
                      double alpha,
                      double beta,
                      RNG& rng) {
      static const std::string function = "beta_binomial_rng";

      check_nonnegative(function, "Population size parameter", N);
      check_positive_finite(function,
                            "First prior sample size parameter", alpha);
      check_positive_finite(function,
                            "Second prior sample size parameter", beta);

      double p = beta_rng(alpha, beta, rng);
      return binomial_rng(N, p, rng);
    }


  }
}
#endif
