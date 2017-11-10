#ifndef STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_RNG_HPP

#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/inv_logit.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>

namespace stan {
  namespace math {

    /**
     * Return pseudorandom Bernoulli draw with specified chance of success 
     * using the specified random number generator.
     *
     * @tparam RNG type of random number generator
     * @param theta chance of success parameter
     * @param rng random number generator
     * @return Bernoulli random variate 
     * @throw std::domain_error if probability parameter is invalid.
     */
    template <class RNG>
    inline int
    bernoulli_rng(double theta,
                  RNG& rng) {
      using boost::variate_generator;
      using boost::bernoulli_distribution;

      static const char* function("bernoulli_rng");

      check_finite(function, "Probability parameter", theta);
      check_bounded(function, "Probability parameter", theta, 0, 1);

      variate_generator<RNG&, bernoulli_distribution<> >
        bernoulli_rng(rng, bernoulli_distribution<>(theta));
      return bernoulli_rng();
    }

  }
}
#endif
