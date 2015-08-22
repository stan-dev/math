#ifndef STAN_MATH_PRIM_SCAL_FUN_INV_PHI_HPP
#define STAN_MATH_PRIM_SCAL_FUN_INV_PHI_HPP

#include <boost/math/tools/promotion.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <stan/math/prim/scal/fun/constants.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_bounded.hpp>
#include <stan/math/prim/scal/fun/Phi.hpp>
#include <stan/math/prim/scal/fun/log1m.hpp>

namespace stan {
  namespace math {

    /**
     * The inverse of the unit normal cumulative distribution function.
     *
     * The return value for a specified input probability, $p$, is
     * the unit normal variate, $x$, such that
     *
     * \f$\Phi(x) = \int_{-\infty}^x \mbox{\sf Norm}(x|0, 1) \ dx = p\f$
     *
     * Algorithm first derived in 2003 by Peter Jon Aklam at
     * http://home.online.no/~pjacklam/notes/invnorm/
     *
     * @param p Argument between 0 and 1.
     * @return Real number
     */
    template <typename T>
    inline typename boost::math::tools::promote_args<T>::type
    inv_Phi(const T p) {
      // overridden in fvar and var, so can hard-code
      // here for scalars only
      using stan::math::check_not_nan;

      check_not_nan("inv_Phi",  "p", p);
      stan::math::check_bounded<T, double, double>
        ("inv_Phi", "Probability variable", p, 0, 1);

      // if T is not double, the gradient will overflow before inv_Phi does
      if (p < 8e-311)
        return NEGATIVE_INFTY;
      if (p == 1)
        return INFTY;

      const double a[6] = {
        -3.969683028665376e+01,  2.209460984245205e+02,
        -2.759285104469687e+02,  1.383577518672690e+02,
        -3.066479806614716e+01,  2.506628277459239e+00
      };
      const double b[5] = {
        -5.447609879822406e+01,  1.615858368580409e+02,
        -1.556989798598866e+02,  6.680131188771972e+01,
        -1.328068155288572e+01
      };
      const double c[6] = {
        -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00,
        4.374664141464968e+00,  2.938163982698783e+00
      };
      const double d[4] = {
        7.784695709041462e-03,  3.224671290700398e-01,
        2.445134137142996e+00,  3.754408661907416e+00
      };

      const double p_low = 0.02425;
      const double p_high = 0.97575;
      T q, r, x;

      if ( (p_low <= p) && (p <= p_high) ) {  // central region
        q = p - 0.5;
        r = q * q;
        x = (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
            (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1.0);
      } else if (p < p_low) {  // lower region
        q = std::sqrt(-2.0*std::log(p));
        x = (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
             ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
      } else {  // upper region
        q = std::sqrt(-2.0*stan::math::log1m(p));
        x = -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
              ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1.0);
      }

      if (x < 37.6) {  // gradient blows up past here anyway
        T e, u;
        e = stan::math::Phi(x) - p;
        u = e * SQRT_2_TIMES_SQRT_PI * std::exp(0.5 * x * x);
        x = x - u / (1.0 + 0.5 * x * u);
      }

      return x;
    }

  }
}

#endif
