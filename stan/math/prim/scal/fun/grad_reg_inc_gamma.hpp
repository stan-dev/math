#ifndef STAN_MATH_PRIM_SCAL_FUN_GRAD_REG_INC_GAMMA_HPP
#define STAN_MATH_PRIM_SCAL_FUN_GRAD_REG_INC_GAMMA_HPP

#include <stan/math/prim/scal/fun/gamma_p.hpp>
#include <stan/math/prim/scal/fun/gamma_q.hpp>
#include <cmath>
#include <stdexcept>

namespace stan {
  namespace math {

    // Gradient of the regularized incomplete gamma functions igamma(a, g)
    // Precomputed values
    // g   = boost::math::tgamma(a)
    // dig = boost::math::digamma(a)
    template<typename T>
    T grad_reg_inc_gamma(T a, T z, T g, T dig, T precision = 1e-6) {
      using boost::math::isinf;
      using std::domain_error;
      using std::exp;
      using std::fabs;
      using std::log;

      T l = log(z);
      if (z >= a && z >= (T)8) {
        // large values of z, compute gradient based on the asymptotic expansion
        // http://dlmf.nist.gov/8.11#E2
        T S = 0;
        T fac = a-1;  // falling_factorial(a-1, k)
        T dfac = 1;   // d/da[falling_factorial(a-1, k)]
        T zpow = z;   // z ** k
        T delta = dfac / zpow;

        for (int k = 1; k < 10; ++k) {
          S += delta;

          zpow *= z;
          dfac = (a-1-k)*dfac + fac;
          fac *= a-1-k;
          delta = dfac / zpow;

          if (isinf(delta))
            throw domain_error("gradRegIncGamma not converging");
        }

        return gamma_q(a, z) * (l - dig) + exp(-z + (a-1)*l) * S / g;
      } else {
        // gradient of series expansion http://dlmf.nist.gov/8.7#E3

        T S = 0;
        T s = 1;
        int k = 0;
        T delta = s / (a * a);
        while (fabs(delta) > precision) {
          S += delta;
          ++k;
          s *= - z / k;
          delta = s / ((k + a) * (k + a));
          if (isinf(delta))
            throw domain_error("gradRegIncGamma not converging");
        }
        return gamma_p(a, z) * ( dig - l ) + exp( a * l ) * S / g;
      }
    }

  }
}
#endif
