#ifndef STAN_MATH_PRIM_SCAL_PROB_UNIFORM_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_UNIFORM_RNG_HPP

#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_greater.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <string>

namespace stan {
  namespace math {

    template <class RNG>
    inline double
    uniform_rng(double alpha,
                double beta,
                RNG& rng) {
      using boost::variate_generator;
      using boost::random::uniform_real_distribution;

      static const std::string function = "uniform_rng";

      check_finite(function, "Lower bound parameter", alpha);
      check_finite(function, "Upper bound parameter", beta);
      check_greater(function, "Upper bound parameter", beta, alpha);

      variate_generator<RNG&, uniform_real_distribution<> >
        uniform_rng(rng, uniform_real_distribution<>(alpha, beta));
      return uniform_rng();
    }

  }
}
#endif
