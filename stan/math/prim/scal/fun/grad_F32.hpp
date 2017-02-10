#ifndef STAN_MATH_PRIM_SCAL_FUN_GRAD_F32_HPP
#define STAN_MATH_PRIM_SCAL_FUN_GRAD_F32_HPP

#include <cmath>

namespace stan {
  namespace math {

    /**
     * Gradients of the hypergeometric function, 3F2.
     *
     * The generalized hypergeometric function is a power series. This
     * implementation computes the gradients using derivatives of the
     * power series directly and stopping when
     * the series converges to within <code>precision</code> or takes
     * <code>max_steps</code>.
     *
     * Although some convergence conditions and divergent conditions are known,
     * this function does not check the inputs for known convergent conditions.
     * Some convergence conditions are listed with the Hypergeometric
     * function itself (F32).
     *
     * @tparam T type of arguments and result
     * @param[out] g g pointer to array of six values of type T, result.
     * @param[in] a1 a1 see generalized hypergeometric function definition.
     * @param[in] a2 a2 see generalized hypergeometric function definition.
     * @param[in] a3 a3 see generalized hypergeometric function definition. 
     * @param[in] b1 b1 see generalized hypergeometric function definition. 
     * @param[in] b2 b2 see generalized hypergeometric function definition. 
     * @param[in] z z see generalized hypergeometric function definition. 
     * @param[in] precision precision of the infinite sum. defaults to 1e-6
     * @param[in] max_steps number of steps to take. defaults to 10000
     */
    template<typename T>
    void grad_F32(T* g, const T& a1, const T& a2, const T& a3, const T& b1, 
        const T& b2, const T& z, const T& precision = 1e-6, int max_steps=1e5) {
      using std::log;
      using std::fabs;
      using std::exp;
      using std::isnan;

      T gOld[6];

      for (T *q = g; q != g + 6; ++q) *q = 0.0;
      for (T *q = gOld; q != gOld + 6; ++q) *q = 0.0;

      T tOld = 1.0;
      T tNew = 0.0;

      T logT = 0.0;

      T logZ = log(z);

      int k = 0;
      T p = 0.0;
      do {
        p = (a1 + k) / (b1 + k) * (a2 + k) / (b2 + k) * (a3 + k) / (1 + k);

        if (isnan(p) || p == 0) 
          break;

        logT += log(fabs(p)) + logZ;
        tNew = sign(p) * exp(logT);

        gOld[0] = tNew * (gOld[0] / tOld + 1.0 / (a1 + k));
        gOld[1] = tNew * (gOld[1] / tOld + 1.0 / (a2 + k));
        gOld[2] = tNew * (gOld[2] / tOld + 1.0 / (a3 + k));

        gOld[3] = tNew * (gOld[3] / tOld - 1.0 / (b1 + k));
        gOld[4] = tNew * (gOld[4] / tOld - 1.0 / (b2 + k));

        gOld[5] = tNew * (gOld[5] / tOld + 1.0 / z);

        for (int i = 0; i < 6; ++i) g[i] += gOld[i];

        tOld = tNew;

        ++k;
        if (k >= max_steps) {
          domain_error("grad_F32", "k (internal counter)", max_steps, 
            "exceeded ", " iterations, hypergeometric function gradient " 
            "did not converge.");
        }
      } while (tNew > precision);
    }

  }
}
#endif
