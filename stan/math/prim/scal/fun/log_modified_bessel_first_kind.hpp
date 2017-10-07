#ifndef STAN_MATH_PRIM_SCAL_FUN_LOG_MODIFIED_BESSEL_FIRST_KIND_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LOG_MODIFIED_BESSEL_FIRST_KIND_HPP

#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_nonnegative.hpp>
#include <stan/math/prim/scal/err/check_greater_or_equal.hpp>
#include <stan/math/prim/scal/fun/is_inf.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <stan/math/prim/scal/fun/lgamma.hpp>
#include <stan/math/prim/scal/fun/log_sum_exp.hpp>
#include <stan/math/prim/scal/fun/log1p.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <limits>

namespace stan {
  namespace math {

    /* Log of the modified Bessel function of the first kind,
     * which is better known as the Bessel I function. See
     * modified_bessel_first_kind.hpp for the function definition.
     *
     * @param v Order, can be a non-integer but must be at least -1
     * @param z Real non-negative number
     * @throws if either v or z is NaN
     * @throws if z is negative
     * @throws if v is less than -1
     * @return log of Bessel I function
     */
    template <typename T1, typename T2>
    inline typename boost::math::tools::promote_args<T1, T2, double>::type
    log_modified_bessel_first_kind(const T1 v, const T2 z) {
      check_not_nan("log_modified_bessel_first_kind", "v", v);
      check_not_nan("log_modified_bessel_first_kind", "z", z);
      check_nonnegative("log_modified_bessel_first_kind", "z", z);
      check_greater_or_equal("log_modified_bessel_first_kind", "v", v, -1);

      using std::log;
      using stan::math::log_sum_exp;
      using stan::math::log1p;
      using stan::math::lgamma;
      using stan::math::square;

      typedef typename boost::math::tools::promote_args<T1, T2, double>::type T;

      if (z == 0) {
       if (v == 0) return 0.0;
       if (v > 0) return -std::numeric_limits<T>::infinity();
       return std::numeric_limits<T>::quiet_NaN();
      }
      if (stan::math::is_inf(z)) return z;
      if (z > 100) {
        // Boost does something like this in asymptotic_bessel_i_large_x
        T lim = (4 * square(v) + 10) / (8 * z);
        lim *= lim;
        lim *= lim;
        lim /= 24;
        if (lim < (std::numeric_limits<double>::epsilon() * 10)) {
          T s = 1;
          T mu = 4 * square(v);
          T ex = 8 * z;
          T num = mu - 1;
          T denom = ex;
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

      T log_half_z = log(0.5 * z);
      T lgam = lgamma(v + 1.0);
      T lcons = (2.0 + v) * log_half_z;
      T out;
      if (v > -1) {
        out = log_sum_exp(v * log_half_z - lgam, lcons - lgamma(v + 2));
        lgam += log(v + 1.0);
      } else {  // lgam is currently NaN
        out = lcons;
        lgam = 0;
      }

      double m = 2;
      T lfac = 0;
      T old_out;
      do {
        lfac += log(m);
        lgam += log(v + m);
        lcons += 2 * log_half_z;
        old_out = out;
        out = log_sum_exp(out, lcons - lfac - lgam);  // underflows eventually
        m++;
      } while (out > old_out || out < old_out);
      return out;
    }

  }
}

#endif
