#ifndef STAN_MATH_PRIM_SCAL_PROB_CHI_SQUARE_RNG_HPP
#define STAN_MATH_PRIM_SCAL_PROB_CHI_SQUARE_RNG_HPP

#include <boost/random/chi_squared_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stan/math/prim/scal/err/check_consistent_sizes.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive_finite.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/fun/multiply_log.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <stan/math/prim/scal/fun/gamma_p.hpp>
#include <stan/math/prim/scal/fun/digamma.hpp>
#include <stan/math/prim/scal/fun/grad_reg_inc_gamma.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>

namespace stan {
  namespace math {

    /**
     * Return a pseudorandom chi squared variate with the nu degrees of freedom
     * using the specified random number generator.
     *
     * @tparam RNG class of random number generator
     * @param nu positive degrees of freedom parameter
     * @param rng random number generator
     * @return chi squared random variate
     * @throw std::domain_error if nu is nonpositive
     */
    template <class RNG>
    inline double
    chi_square_rng(double nu,
                   RNG& rng) {
      using boost::variate_generator;
      using boost::random::chi_squared_distribution;

      static const char* function("chi_square_rng");

      check_positive_finite(function, "Degrees of freedom parameter", nu);

      variate_generator<RNG&, chi_squared_distribution<> >
        chi_square_rng(rng, chi_squared_distribution<>(nu));
      return chi_square_rng();
    }

  }
}
#endif
