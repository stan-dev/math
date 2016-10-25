#ifndef STAN_MATH_PRIM_SCAL_PROB_BETA_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_BETA_RNG_HPP

#include <boost/math/special_functions/gamma.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_less_or_equal.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/lbeta.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#include <stan/math/prim/scal/fun/grad_reg_inc_beta.hpp>
#include <stan/math/prim/scal/fun/inc_beta.hpp>
#include <stan/math/prim/mat/prob/dirichlet_rng.hpp>

namespace stan {
  namespace math {

    template <class RNG>
    inline double
    beta_rng(double alpha,
             double beta,
             RNG& rng) {
      using Eigen::Dynamic;
      using Eigen::Matrix;
      using stan::math::dirichlet_rng;
      
      static const char* function("beta_rng");

      check_positive_finite(function, "shape parameter", alpha);
      check_positive_finite(function, "shape parameter", beta);

      Matrix<double, Dynamic, 1> dirichlet_params(2, 1);
      dirichlet_params << alpha, beta;
      return dirichlet_rng(dirichlet_params, rng)[0];
    }

  }
}
#endif
