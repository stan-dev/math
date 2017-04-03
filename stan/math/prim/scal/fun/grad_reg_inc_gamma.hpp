#ifndef STAN_MATH_PRIM_SCAL_FUN_GRAD_REG_INC_GAMMA_HPP
#define STAN_MATH_PRIM_SCAL_FUN_GRAD_REG_INC_GAMMA_HPP

#include <stan/math/prim/scal/err/domain_error.hpp>
#include <stan/math/prim/scal/fun/gamma_p.hpp>
#include <stan/math/prim/scal/fun/gamma_q.hpp>
#include <stan/math/prim/scal/fun/is_inf.hpp>
#include <stan/math/prim/scal/fun/square.hpp>
#include <cmath>

namespace stan {
  namespace math {

    /**
     * Gradient of the regularized incomplete gamma functions igamma(a, z)
     *
     * For small z, the gradient is computed via the series expansion;
     * for large z, the series is numerically inaccurate due to cancellation
     * and the asymptotic expansion is used.
     *
     * @param a   shape parameter, a > 0
     * @param z   location z >= 0
     * @param g   boost::math::tgamma(a) (precomputed value)
     * @param dig boost::math::digamma(a) (precomputed value)
     * @param precision required precision; applies to series expansion only
     * @param max_steps number of steps to take. defaults to 10000
     * @throw throws std::domain_error if not converged after max_steps
     *   or increment overflows to inf.
     *
     * For the asymptotic expansion, the gradient is given by:
       \f[
       \begin{array}{rcl}
       \Gamma(a, z) & = & z^{a-1}e^{-z} \sum_{k=0}^N \frac{(a-1)_k}{z^k} \qquad , z \gg a\\
       Q(a, z) & = & \frac{z^{a-1}e^{-z}}{\Gamma(a)} \sum_{k=0}^N \frac{(a-1)_k}{z^k}\\
       (a)_k & = & (a)_{k-1}(a-k)\\
       \frac{d}{da} (a)_k & = & (a)_{k-1} + (a-k)\frac{d}{da} (a)_{k-1}\\
       \frac{d}{da}Q(a, z) & = & (log(z) - \psi(a)) Q(a, z)\\
       && + \frac{z^{a-1}e^{-z}}{\Gamma(a)} \sum_{k=0}^N \left(\frac{d}{da} (a-1)_k\right) \frac{1}{z^k}
       \end{array}
       \f]
     */
    template<typename T>
    T grad_reg_inc_gamma(T a, T z, T g, T dig, double precision = 1e-6,
        int max_steps = 1e5) {
      using std::domain_error;
      using std::exp;
      using std::fabs;
      using std::log;

      T l = log(z);
      if (z >= a && z >= 8) {
        // asymptotic expansion http://dlmf.nist.gov/8.11#E2
        T S = 0;
        T a_minus_one_minus_k = a - 1;
        T fac = a_minus_one_minus_k;  // falling_factorial(a-1, k)
        T dfac = 1;   // d/da[falling_factorial(a-1, k)]
        T zpow = z;   // z ** k
        T delta = dfac / zpow;

        for (int k = 1; k < 10; ++k) {
          a_minus_one_minus_k -= 1;

          S += delta;

          zpow *= z;
          dfac = a_minus_one_minus_k * dfac + fac;
          fac *= a_minus_one_minus_k;
          delta = dfac / zpow;

          if (is_inf(delta))
            stan::math::domain_error("grad_reg_inc_gamma",
                                     "is not converging", "", "");
        }

        return gamma_q(a, z) * (l - dig) + exp(-z + (a - 1) * l) * S / g;
      } else {
        // gradient of series expansion http://dlmf.nist.gov/8.7#E3

        T S = 0;
        T s = 1;
        int k = 0;
        T delta = s / square(a);
        while (fabs(delta) > precision) {
          S += delta;
          ++k;
          s *= - z / k;
          delta = s / square(k + a);
          if (k >= max_steps)
            stan::math::domain_error("grad_reg_inc_gamma",
              "k (internal counter)",
              max_steps, "exceeded ",
              " iterations, gamma function gradient did not converge.");
          if (is_inf(delta))
            stan::math::domain_error("grad_reg_inc_gamma",
                                     "is not converging", "", "");
        }
        return gamma_p(a, z) * ( dig - l ) + exp( a * l ) * S / g;
      }
    }

  }
}
#endif
