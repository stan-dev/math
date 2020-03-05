#ifndef STAN_MATH_PRIM_FUN_INC_BETA_DDB_HPP
#define STAN_MATH_PRIM_FUN_INC_BETA_DDB_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/fabs.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/inc_beta_dda.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T>
T inc_beta_dda(T a, T b, T z, T digamma_a, T digamma_ab);

/**
 * Returns the partial derivative of the regularized
 * incomplete beta function, I_{z}(a, b) with respect to b.
 * The power series used to compute the derivative tends to
 * converge slowly when a and b are large, especially if z
 * approaches 1.  The implementation will throw an exception
 * if the series have not converged within 100,000 iterations.
 * The current implementation has been tested for values
 * of a and b up to 12500 and z = 0.999.
 *
 * @tparam T scalar types of arguments
 * @param a first argument
 * @param b second argument
 * @param z upper bound of the integral
 * @param digamma_b value of digamma(b)
 * @param digamma_ab value of digamma(b)
 * @return partial derivative of the incomplete beta with respect to b
 *
 * @pre a >= 0
 * @pre b >= 0
 * @pre 0 <= z <= 1
 */
template <typename T>
T inc_beta_ddb(T a, T b, T z, T digamma_b, T digamma_ab) {
  using std::fabs;
  using std::log;
  using std::pow;

  if (b > a) {
    if ((0.1 < z && z <= 0.75 && b > 500) || (0.01 < z && z <= 0.1 && b > 2500)
        || (0.001 < z && z <= 0.01 && b > 1e5)) {
      return -inc_beta_dda(b, a, 1 - z, digamma_b, digamma_ab);
    }
  }

  if ((z > 0.75 && a < 500) || (z > 0.9 && a < 2500) || (z > 0.99 && a < 1e5)
      || (z > 0.999)) {
    return -inc_beta_dda(b, a, 1 - z, digamma_b, digamma_ab);
  }

  double threshold = 1e-10;

  const T a_plus_b = a + b;
  const T a_plus_1 = a + 1;

  // Common prefactor to regularize numerator and denominator
  T prefactor = pow(a_plus_1 / a_plus_b, 3);

  T sum_numer = digamma_ab * prefactor;
  T sum_denom = prefactor;

  T summand = prefactor * z * a_plus_b / a_plus_1;

  T k = 1;
  digamma_ab += inv(a_plus_b);

  while (fabs(summand) > threshold) {
    sum_numer += digamma_ab * summand;
    sum_denom += summand;

    summand *= (1 + (a_plus_b) / k) * (1 + k) / (1 + a_plus_1 / k);
    digamma_ab += inv(a_plus_b + k);
    ++k;
    summand *= z / k;

    if (k > 1e5) {
      throw_domain_error("inc_beta_ddb",
                         "did not converge within 100000 iterations", "", "");
    }
  }

  return inc_beta(a, b, z) * (log(1 - z) - digamma_b + sum_numer / sum_denom);
}

}  // namespace math
}  // namespace stan
#endif
