#ifndef STAN_MATH_PRIM_FUN_INC_BETA_DDZ_HPP
#define STAN_MATH_PRIM_FUN_INC_BETA_DDZ_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/log1m.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the partial derivative of the regularized
 * incomplete beta function, I_{z}(a, b) with respect to z.
 *
 * @tparam T scalar types of arguments
 * @param a first argument
 * @param b second argument
 * @param z upper bound of the integral
 * @return partial derivative of the incomplete beta with respect to z
 *
 * @pre a > 0
 * @pre b > 0
 * @pre 0 < z <= 1
 */
template <typename T>
T inc_beta_ddz(T a, T b, T z) {
  using std::exp;
  using std::log;
  return exp((b - 1) * log1m(z) + (a - 1) * log(z) + lgamma(a + b) - lgamma(a)
             - lgamma(b));
}

template <>
inline double inc_beta_ddz(double a, double b, double z) {
  using boost::math::ibeta_derivative;
  return ibeta_derivative(a, b, z);
}

}  // namespace math
}  // namespace stan
#endif
