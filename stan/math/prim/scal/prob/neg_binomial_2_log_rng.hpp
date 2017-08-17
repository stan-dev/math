#ifndef STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_2_LOG_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_NEG_BINOMIAL_2_LOG_RNG_HPP

#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_less.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <string>

namespace stan {
  namespace math {

    template <class RNG>
    inline int
    neg_binomial_2_log_rng(double eta,
                           double phi,
                           RNG& rng) {
      using boost::variate_generator;
      using boost::random::poisson_distribution;
      using boost::random::gamma_distribution;

      static const std::string function = "neg_binomial_2_log_rng";

      check_finite(function, "Log-location parameter", eta);
      check_positive_finite(function, "Precision parameter", phi);

      double exp_eta_div_phi = std::exp(eta)/phi;

      // gamma_rng params must be positive and finite
      check_positive_finite(function,
                            "Exponential of the log-location parameter "
                            "divided by the precision parameter",
                            exp_eta_div_phi);

      double rng_from_gamma =
        variate_generator<RNG&, gamma_distribution<> >
        (rng, gamma_distribution<>(phi, exp_eta_div_phi))();

      // same as the constraints for poisson_rng
      check_less(function,
                 "Random number that came from gamma distribution",
                 rng_from_gamma, POISSON_MAX_RATE);
      check_not_nan(function,
                    "Random number that came from gamma distribution",
                    rng_from_gamma);
      check_nonnegative(function,
                        "Random number that came from gamma distribution",
                        rng_from_gamma);

      return variate_generator<RNG&, poisson_distribution<> >
        (rng, poisson_distribution<>(rng_from_gamma))();
    }

  }
}
#endif
