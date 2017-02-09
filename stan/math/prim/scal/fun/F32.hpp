#ifndef STAN_MATH_PRIM_SCAL_FUN_F32_HPP
#define STAN_MATH_PRIM_SCAL_FUN_F32_HPP

#include <stan/math/prim/scal/fun/sign.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Generalized hypergeometric function, 3F2.
     *
     * The generalized hypergeometric function is a power series. This
     * implementation computes the power series directly stopping when
     * the series converges to within <code>precision</code> or takes
     * <code>max_steps</code>.
     *
     * Although some convergence conditions and divergent conditions are known,
     * this function does not check the inputs for known converent conditions.
     *
     * This function will converge if:
     *   - <code>a1</code>, <code>a2</code>, or <code>a3</code> is a
     *     non-positive integer
     *   - <code>b1</code> or <code>b2</code> is a non-positive integer
     *   - <code>z</code> is less than 1
     * This function will diverge if <code>z</code> is greater than 1.
     * When <code>z</code> is 1, which is the case for the beta binomial
     * cdf, it is hard to determine.
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
          const T& z, double precision = 1e-6, int max_steps = 10000) {
      using std::exp;
      using std::log;
      using std::fabs;

      T F = 1.0;
      T tNew = 0.0;
      T logT = 0.0;
      T logZ = log(z);

      int k = 0;
      T p = 0;
      do {
        p = (a1 + k) * (a2 + k) * (a3 + k) / ((b1 + k) * (b2 + k) * (k + 1));

        if (p == 0)
          break;

        logT += sign(p) * log(fabs(p)) + logZ;
        tNew = exp(logT);
        F += tNew;
        ++k;
      } while (tNew > precision && k < max_steps);
      return F;
    }

  }
}
#endif
