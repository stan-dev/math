#ifndef STAN_MATH_PRIM_SCAL_PROB_BETA_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BETA_RNG_HPP

#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/log_sum_exp.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/lbeta.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/scal/fun/inc_beta.hpp>

namespace stan {
  namespace math {

    /**
     * Return a pseudorandom Beta variate with the supplied success and failure
     * parameters and specified random number generator.
     *
     * @tparam RNG class of random number generator
     * @param alpha positive finite success parameter
     * @param beta positive finite failure parameter
     * @param rng random number generator
     * @return Beta random variate
     * @throw std::domain_error if alpha or beta is nonpositive
     */
    template <class RNG>
    inline double
    beta_rng(double alpha,
             double beta,
             RNG& rng) {
      using boost::variate_generator;
      using boost::random::gamma_distribution;
      using boost::random::uniform_real_distribution;
      using std::log;
      using std::exp;
      static const char* function("beta_rng");
      check_positive_finite(function, "First shape parameter", alpha);
      check_positive_finite(function, "Second shape parameter", beta);

      // If alpha and beta are large, trust the usual ratio of gammas
      // method for generating beta random variables. If any parameter
      // is small, work in log space and use Marsaglia and Tsang's trick
      if (alpha > 1.0 && beta > 1.0) {
        variate_generator<RNG&, gamma_distribution<> >
          rng_gamma_alpha(rng, gamma_distribution<>(alpha, 1.0));
        variate_generator<RNG&, gamma_distribution<> >
          rng_gamma_beta(rng, gamma_distribution<>(beta, 1.0));
        double a = rng_gamma_alpha();
        double b = rng_gamma_beta();
        return a / (a + b);
      } else {
        variate_generator<RNG&, uniform_real_distribution<> >
          uniform_rng(rng, uniform_real_distribution<>(0.0, 1.0));
        variate_generator<RNG&, gamma_distribution<> >
          rng_gamma_alpha(rng, gamma_distribution<>(alpha + 1, 1.0));
        variate_generator<RNG&, gamma_distribution<> >
          rng_gamma_beta(rng, gamma_distribution<>(beta + 1, 1.0));
        double log_a = log(uniform_rng()) / alpha + log(rng_gamma_alpha());
        double log_b = log(uniform_rng()) / beta + log(rng_gamma_beta());
        double log_sum = log_sum_exp(log_a, log_b);
        return exp(log_a - log_sum);
      }
    }

  }
}
#endif
