#ifndef STAN_MATH_PRIM_SCAL_FUN_F32_HPP
#define STAN_MATH_PRIM_SCAL_FUN_F32_HPP

#include <stan/math/prim/scal/fun/sign.hpp>
#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/scal/fun/is_nan.hpp>
#include <stan/math/prim/scal/fun/is_inf.hpp>
#include <stan/math/prim/scal/err/check_3F2_converges.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Hypergeometric function, 3F2.
     *
     * Calculate the hypergeometric function (3F2) as the power series
     * directly to within <code>precision</code> or until
     * <code>max_steps</code> terms.
     *
     * This function does not have a closed form but will converge if:
     *   - <code>|z|</code> is less than 1
     *   - <code>|z|</code> is equal to one and <code>b1 + b2 < a1 + a2 + a3</code>
     * This function is a rational polynomial if
     *   - <code>a1</code>, <code>a2</code>, or <code>a3</code> is a
     *     non-positive integer
     * This function can be treated as a rational polynomial if
     *   - <code>b1</code> or <code>b2</code> is a non-positive integer
     *     and the series is terminated prior to the final term.
     *
     * @tparam T type of arguments and result
     * @param[in] a1 a1 (always called with 1 from beta binomial cdfs)
     * @param[in] a2 a2 (always called with a2 > 1)
     * @param[in] a3 a3 (always called with int a3 <= 0)
     * @param[in] b1 b1 (always called with int b1 < |a3|)
     * @param[in] b2 b2 (always <= 1)
     * @param[in] z z (is always called with 1 from beta binomial cdfs)
     * @param[in] precision precision of the infinite sum. defaults to 1e-6
     * @param[in] max_steps number of steps to take. defaults to 10000
     */
    template<typename T>
    T F32(const T& a1, const T& a2, const T& a3, const T& b1, const T& b2,
          const T& z, double precision = 1e-6, int max_steps = 1e5) {

      check_3F2_converges("F32", a1, a2, a3, b1, b2, z);

      using std::exp;
      using std::log;
      using std::fabs;
      using stan::math::is_nan;

      T F = 1.0;
      T tNew = 0.0;
      T logT = 0.0;
      T logZ = log(z);

      int k = 0;
      bool T_is_negative = false;
      T p = 0.0;
      do {
        p = (a1 + k) * (a2 + k) * (a3 + k) / ((b1 + k) * (b2 + k) * (k + 1));
        if (is_nan(p) || p == 0)
          break;

        logT += log(fabs(p)) + logZ;
        if (p < 0 && T_is_negative) {
          T_is_negative = false;
        } else if (p < 0 && !T_is_negative) {
          T_is_negative = true;
        }
        if (T_is_negative)
          tNew = -1 * exp(logT);
        else
          tNew = exp(logT);
        F += tNew;

        ++k;
        if (k >= max_steps) {
          domain_error("F32", "k (internal counter)", max_steps, "exceeded ",
            " iterations, hypergeometric function did not converge.");
        }
        if (is_inf(F)) {
          domain_error("F32", "F (output)", F,
            "overflow ", " hypergeometric function did not converge.");
        }
      } while (fabs(tNew) > precision);
      return F;
    }

  }
}
#endif
