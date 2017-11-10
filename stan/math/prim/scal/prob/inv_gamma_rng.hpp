#ifndef STAN_MATH_PRIM_SCAL_PROB_INV_GAMMA_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_INV_GAMMA_RNG_HPP

#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/gamma_q.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/meta/length.hpp>
#include <stan/math/prim/scal/meta/VectorBuilder.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
  namespace math {

    template <class RNG>
    inline double
    inv_gamma_rng(double alpha,
                  double beta,
                  RNG& rng) {
      using boost::variate_generator;
      using boost::random::gamma_distribution;

      static const char* function("inv_gamma_rng");

      check_positive_finite(function, "Shape parameter", alpha);
      check_positive_finite(function, "Scale parameter", beta);

      variate_generator<RNG&, gamma_distribution<> >
        gamma_rng(rng, gamma_distribution<>(alpha, 1 / beta));
      return 1 / gamma_rng();
    }

  }
}
#endif
