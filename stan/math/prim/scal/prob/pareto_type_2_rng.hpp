#ifndef STAN_MATH_PRIM_SCAL_PROB_PARETO_TYPE_2_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_PARETO_TYPE_2_RNG_HPP

#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/scal/prob/uniform_rng.hpp>
#include <boost/random/variate_generator.hpp>
#include <cmath>
#include <string>

namespace stan {
  namespace math {

    template <class RNG>
    inline double
    pareto_type_2_rng(double mu,
                      double lambda,
                      double alpha,
                      RNG& rng) {
      static const std::string function = "pareto_type_2_rng";

      check_positive(function, "scale parameter", lambda);
      double uniform_01 = uniform_rng(0.0, 1.0, rng);
      return (std::pow(1.0 - uniform_01, -1.0 / alpha) - 1.0) * lambda + mu;
    }

  }
}
#endif
