#ifndef STAN_MATH_PRIM_SCAL_FUN_F32_HPP
#define STAN_MATH_PRIM_SCAL_FUN_F32_HPP

#include <stan/math/prim/scal/fun/sign.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Generalized hypergeometric function, 3F2.
     * Implementation isn't general and doesn't work for all values
     * of inputs.
     *
     * @tparam T type of arguments and result
     * @param a1 a1 (always called with 1 from beta binomial cdfs)
     * @param a2 a2 (always called with a2 > 1)
     * @param a3 a3 (always called with int a3 <= 0)
     * @param b1 b1 (always called with int b1 < |a3|)
     * @param b2 b2 (always <= 1)
     * @param z z (is always called with 1 from beta binomial cdfs)
     * @param precision precision of the infinite sum. defaults to 1e-6
     * @param max_steps number of steps to take. defaults to 10000
     */
    template<typename T>
    T F32(T a1, T a2, T a3,
          T b1, T b2,
          T z,
          double precision = 1e-6,
          int max_steps = 10000) {
      using std::exp;
      using std::log;
      using std::fabs;

      T F = 1.0;
      T tNew = 0.0;
      T logT = 0.0;
      T logZ = log(z);

      int k = 0;
      do {
        T p = (a1 + k) * (a2 + k) * (a3 + k) / ((b1 + k) * (b2 + k) * (k + 1));

        // If a1, a2, or a3 is a negative integer then the series terminates
        // after a finite number of interations
        if (p == 0)
          break;

        logT += sign(p) * log(fabs(p)) + logZ;
        tNew = exp(logT);
        F += tNew;
        ++k;
      } while (fabs(tNew) > precision && k < max_steps);
      return F;
    }

  }
}
#endif
