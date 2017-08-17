#ifndef STAN_MATH_PRIM_SCAL_PROB_PARETO_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_PARETO_RNG_HPP

#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>
#include <string>

namespace stan {
  namespace math {

    template <class RNG>
    inline double
    pareto_rng(double y_min,
               double alpha,
               RNG& rng) {
      using boost::variate_generator;
      using boost::exponential_distribution;

      static const std::string function = "pareto_rng";

      check_positive_finite(function, "Scale parameter", y_min);
      check_positive_finite(function, "Shape parameter", alpha);

      variate_generator<RNG&, exponential_distribution<> >
        exp_rng(rng, exponential_distribution<>(alpha));
      return y_min * std::exp(exp_rng());
    }

  }
}
#endif
