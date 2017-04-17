#ifndef STAN_MATH_PRIM_SCAL_PROB_GUMBEL_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_GUMBEL_RNG_HPP

#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
  namespace math {

    /**
     * Return a pseudorandom Gumbel variate with the given location and scale
     * using the specified random number generator.
     *
     * @tparam RNG type of random number generator
     * @param mu location parameter
     * @param beta positive scale parameter
     * @param rng random number generator
     * @return Gumbel random variate
     * @throw std::domain_error if mu is infinite or beta is nonpositive.
     */
    template <class RNG>
    inline double
    gumbel_rng(double mu,
               double beta,
               RNG& rng) {
      using boost::variate_generator;
      using boost::uniform_01;

      static const char* function("gumbel_rng");

      check_finite(function, "Location parameter", mu);
      check_positive(function, "Scale parameter", beta);

      variate_generator<RNG&, uniform_01<> >
        uniform01_rng(rng, uniform_01<>());
      return mu - beta * std::log(-std::log(uniform01_rng()));
    }

  }
}
#endif
