#ifndef STAN_MATH_PRIM_SCAL_PROB_BETA_BINOMIAL_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BETA_BINOMIAL_RNG_HPP

#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/binomial_coefficient_log.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/lbeta.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/prob/binomial_rng.hpp>
#include <stan/math/prim/scal/prob/beta_rng.hpp>
#include <stan/math/prim/scal/fun/F32.hpp>
#include <stan/math/prim/scal/fun/grad_F32.hpp>

namespace stan {
  namespace math {

    /**
     * Paranoid wrapper for beta_rng which enforces good output at the cost
     * of introducing an infinite loop.
     *
     * @param alpha success parameter
     * @param beta failure parameter
     * @param rng random number generator
     * @tparam RNG class of random number generator
     * @return sample from beta distribution
     */
    template <class RNG>
    inline double safe_beta_rng(double alpha, double beta, RNG& rng) {
      double p = beta_rng(alpha, beta, rng);
      while (p < 0 || p > 1) {
        p = beta_rng(alpha, beta, rng);
      }
      return p;
    }

    /**
     * Beta-Binomial random number generator. Equivalent to drawing p from Beta(alpha,beta)
     * and then drawing from a Binomial(N,p) distribution.
     *
     * @param N nonnegative population size.
     * @param alpha nonnegative sucess parameter.
     * @param beta nonnegative failure parameter.
     * @param rng random number generator.
     *
     * @tparam RNG type of rng.
     *
     * @return draw from beta-binomial distribution with specified parameters.
     *
     * @throw std::domain_error if N, alpha, or beta is negative.
     */
    template <class RNG>
    inline int
    beta_binomial_rng(int N,
                      double alpha,
                      double beta,
                      RNG& rng) {
      static const char* function("beta_binomial_rng");

      check_nonnegative(function, "Population size parameter", N);
      check_positive_finite(function,
                            "First prior sample size parameter", alpha);
      check_positive_finite(function,
                            "Second prior sample size parameter", beta);

      double p = safe_beta_rng(alpha, beta, rng);
      return binomial_rng(N, p, rng);
    }


  }
}
#endif
