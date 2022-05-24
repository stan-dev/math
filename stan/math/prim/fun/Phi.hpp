#ifndef STAN_MATH_PRIM_FUN_PHI_HPP
#define STAN_MATH_PRIM_FUN_PHI_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erf.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <stan/math/prim/fun/Phi.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>

namespace stan {
namespace math {

/**
 * The unit normal cumulative distribution function.
 *
 * The return value for a specified input is the probability that
 * a random unit normal variate is less than or equal to the
 * specified value, defined by
 *
 * \f$\Phi(x) = \int_{-\infty}^x \mbox{\sf Norm}(x|0, 1) \ dx\f$
 *
 * This function can be used to implement the inverse link function
 * for probit regression.
 *
 * Phi will underflow to 0 below -37.5 and overflow to 1 above 8
 *
 * @param x Argument.
 * @return Probability random sample is less than or equal to argument.
 */
inline double Phi(double x) {
  check_not_nan("Phi", "x", x);
  if (x < -37.5) {
    return 0;
  } else if (x < -5.0) {
    return 0.5 * erfc(-INV_SQRT_TWO * x);
  } else if (x > 8.25) {
    return 1;
  } else {
    return 0.5 * (1.0 + erf(INV_SQRT_TWO * x));
  }
}

/**
 * Structure to wrap Phi() so it can be vectorized.
 *
 * @tparam T type of argument
 * @param x argument
 * @return Unit normal CDF of x.
 */
struct Phi_fun {
  template <typename T>
  static inline auto fun(const T& x) {
    return Phi(x);
  }
};

/**
 * Vectorized version of Phi().
 *
 * @tparam T type of container
 * @param x container
 * @return Unit normal CDF of each value in x.
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_not_var_matrix_t<T>* = nullptr>
inline auto Phi(const T& x) {
  return apply_scalar_unary<Phi_fun, T>::apply(x);
}

}  // namespace math
}  // namespace stan

#endif
