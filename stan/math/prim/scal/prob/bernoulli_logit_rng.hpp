#ifndef STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_LOGIT_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BERNOULLI_LOGIT_RNG_HPP

#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/fun/inv_logit.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <string>

namespace stan {
  namespace math {

    /**
     * A Bernoulli random number generator which takes as its argument the
     * often more convenient logit-parametrization.
     *
     * @tparam RNG Random number generator type.
     * @param t logit-transformed probability parameter.
     * @param rng pseudorandom number generator.
     * @return Bernoulli(logit^{-1}(t)) generated random number, either 0 or 1.
     */
    template <class RNG>
    inline int
    bernoulli_logit_rng(double t, RNG& rng) {
      using boost::variate_generator;
      using boost::bernoulli_distribution;
      using stan::math::inv_logit;

      check_finite("bernoulli_logit_rng",
                   "Logit transformed probability parameter", t);

      variate_generator<RNG&, bernoulli_distribution<> >
        bernoulli_rng(rng, bernoulli_distribution<>(inv_logit(t)));
      return bernoulli_rng();
    }

  }
}
#endif
