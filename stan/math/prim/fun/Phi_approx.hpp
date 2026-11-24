#ifndef STAN_MATH_PRIM_FUN_PHI_APPROX_HPP
#define STAN_MATH_PRIM_FUN_PHI_APPROX_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/functor/apply_scalar_unary.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return an approximation of the unit normal CDF.
 *
 * Computes \f$ \mbox{logit}^{-1}(\gamma x^3 + \beta x) \f$ with the
 * coefficients <code>PHI_APPROX_CUBIC_COEF</code> and
 * <code>PHI_APPROX_LINEAR_COEF</code>, the form of Bowling et al. (2009),
 * http://www.jiem.org/index.php/jiem/article/download/60/27
 *
 * This function can be used to implement the inverse link function
 * for probit regression.
 *
 * @param x Argument.
 * @return Probability random sample is less than or equal to argument.
 */
inline double Phi_approx(double x) {
  using std::pow;
  return inv_logit(PHI_APPROX_CUBIC_COEF * pow(x, 3.0)
                   + PHI_APPROX_LINEAR_COEF * x);
}

/**
 * Return an approximation of the unit normal CDF.
 *
 * @param x argument.
 * @return approximate probability random sample is less than or
 * equal to argument.
 */
inline double Phi_approx(int x) { return Phi_approx(static_cast<double>(x)); }

/**
 * Structure to wrap Phi_approx() so it can be vectorized.
 */
struct Phi_approx_fun {
  /**
   * Return the approximate value of the Phi() function applied to
   * the argument.
   *
   * @tparam T type of argument
   * @param x argument
   * @return approximate value of Phi applied to argument
   */
  template <typename T>
  static inline auto fun(T&& x) {
    return Phi_approx(std::forward<T>(x));
  }
};

/**
 * Return the elementwise application of <code>Phi_approx()</code> to
 * specified argument container.  The return type promotes the
 * underlying scalar argument type to double if it is an integer,
 * and otherwise is the argument type.
 *
 * @tparam T type of container
 * @param x container
 * @return elementwise Phi_approx of container elements
 */
template <
    typename T,
    require_all_not_nonscalar_prim_or_rev_kernel_expression_t<T>* = nullptr,
    require_container_t<T>* = nullptr, require_not_var_matrix_t<T>* = nullptr>
inline auto Phi_approx(T&& x) {
  return apply_scalar_unary<Phi_approx_fun, T>::apply(std::forward<T>(x));
}

}  // namespace math
}  // namespace stan

#endif
