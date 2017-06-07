#ifndef STAN_MATH_PRIM_SCAL_FUN_LOG_MODIFIED_BESSEL_FIRST_KIND_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LOG_MODIFIED_BESSEL_FIRST_KIND_HPP

#include <boost/math/special_functions/bessel.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <stan/math/prim/scal/fun/is_inf.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/log_sum_exp.hpp>
#include <stan/math/prim/scal/fun/log1p.hpp>
#include <limits>

namespace stan {
  namespace math {

    /* Log of the modified Bessel function of the first kind,
     * which is better known as the log of the Bessel I function
     *
     * @param v Order, can be a non-integer but must be at least -1
     * @param z Real non-negative number
     * @throws if either v or z is NaN
     * @throws if z is negative
     * @throws if v is less than -1
     * @return log of Bessel I function
     */
    inline double
    log_modified_bessel_first_kind(const double v, const double z) {
      check_not_nan("log_modified_bessel_first_kind", "v", v);
      check_not_nan("log_modified_bessel_first_kind", "z", z);
      check_nonnegative("log_modified_bessel_first_kind", "z", z);
      check_greater_or_equal("log_modified_bessel_first_kind", "v", v, -1);

      using std::log;
      using stan::math::log_sum_exp;
      using stan::math::log1p;
      using stan::math::lgamma;

      if (z == 0) {
       if (v == 0) return 0.0;
       if (v > 0) return -std::numeric_limits<double>::infinity();
       return std::numeric_limits<double>::quiet_NaN();
      }
      if (stan::math::is_inf(z)) return z;
      if (z > 100) {
        // Boost does something like this in asymptotic_bessel_i_large_x
        double lim = (4 * v * v + 10) / (8 * z);
        lim *= lim;
        lim *= lim;
        lim /= 24;
        if (lim < (std::numeric_limits<double>::epsilon() * 10)) {
          double s = 1;
          double mu = 4 * v * v;
          double ex = 8 * z;
          double num = mu - 1;
          double denom = ex;
          s -= num / denom;
          num *= mu - 9;
          denom *= ex * 2;
          s += num / denom;
          num *= mu - 25;
          denom *= ex * 3;
          s -= num / denom;
          s = z - log(std::sqrt(2 * z * M_PI)) + log(s);
          return s;
        }
      }

      double log_half_z = log(0.5 * z);
      double lgam = lgamma(v + 1);
      double lcons = (2 + v) * log_half_z;
      double out;
      if (v > -1) {
        out = log_sum_exp(v * log_half_z - lgam, lcons - lgamma(v + 2));
        lgam += log(v + 1);
      } else {  // lgam is currently NaN
        out = lcons;
        lgam = 0;
      }

      int m = 2;
      double lfac = 0;
      double old_out = std::numeric_limits<double>::quiet_NaN();
      while (out != old_out) {
        lfac += log(m);
        lgam += log(v + m);
        lcons += 2 * log_half_z;
        old_out = out;
        out = log_sum_exp(out, lcons - lfac - lgam);  // underflows eventually
        m++;
      }
      return out;
    }

  }
}

#endif
