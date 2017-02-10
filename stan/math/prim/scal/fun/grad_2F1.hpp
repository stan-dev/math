#ifndef STAN_MATH_PRIM_SCAL_FUN_GRAD_2F1_HPP
#define STAN_MATH_PRIM_SCAL_FUN_GRAD_2F1_HPP

#include <cmath>

namespace stan {
  namespace math {

    // Gradient of the hypergeometric function 2F1(a, b | c | z)
    // with respect to a and c
    /**
     * Gradient of the hypergeometric function, 2F1(a1, a2, b1, z)
     * with respect to a1 and b1 only.
     *
     * The generalized hypergeometric function is a power series. This
     * implementation computes gradient by computing the power series 
     * directly stopping when the series converges to within 
     * <code>precision</code> or takes <code>max_steps</code>.
     *
     * If more than </code>max_steps</code> are taken without 
     * converging, the function will throw a domain_error.
     *
     * @tparam T type of arguments and result
     * @param[out] gradA11 output argument for partial w.r.t. a1
     * @param[out] gradB1 output argument for partial w.r.t. b1
     * @param[in] a1 a1, see generalized hypergeometric function definition.
     * @param[in] a2 a2, see generalized hypergeometric function definition.
     * @param[in] b1 b1, see generalized hypergeometric function definition.
     * @param[in] z z, see generalized hypergeometric function definition.
     * @param[in] precision precision of the infinite sum. defaults to 1e-6
     * @param[in] max_steps number of steps to take. defaults to 10000
     * @throw @throws std::domain_error if not converged after max_steps
     * 
     */
    template<typename T>
    void grad_2F1(T& gradA1, T& gradB1, const T& a1, const T& a2, 
      const T& b1, const T& z, T precision = 1e-6, int max_steps = 1e5) {
      using std::fabs;
      using std::isnan;

      gradA1 = 0;
      gradB1 = 0;

      T gradA1old = 0;
      T gradB1old = 0;

      int k = 0;
      T tDak = 1.0 / (a1 - 1);

      do {
        const T r = ( (a1 + k) / (b1 + k) ) * ( (a2 + k) / (k + 1) ) * z;
        tDak = r * tDak * (a1 + (k - 1)) / (a1 + k);

        if (isnan(r) || r == 0) 
          break;

        gradA1old = r * gradA1old + tDak;
        gradB1old = r * gradB1old - tDak * ((a1 + k) / (b1 + k));

        gradA1 += gradA1old;
        gradB1 += gradB1old;

        ++k;
        if (k >= max_steps) {
          domain_error("grad_2F1", "k (internal counter)", max_steps, 
            "exceeded ", 
            " iterations, hypergeometric function did not converge.");
        }

      } while (fabs(tDak * (a1 + (k - 1)) ) > precision) 
    }

  }
}
#endif
