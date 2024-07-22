#ifndef STAN_MATH_PRIM_FUN_LMGAMMA_HPP
#define STAN_MATH_PRIM_FUN_LMGAMMA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>

namespace stan {
namespace math {

/**
 * Return the natural logarithm of the multivariate gamma function
 * with the specified dimensions and argument.
 *
 * <p>The multivariate gamma function \f$\Gamma_k(x)\f$ for
 * dimensionality \f$k\f$ and argument \f$x\f$ is defined by
 *
 * <p>\f$\Gamma_k(x) = \pi^{k(k-1)/4} \, \prod_{j=1}^k \Gamma(x + (1 - j)/2)\f$,
 *
 * where \f$\Gamma()\f$ is the gamma function.
 *
   \f[
   \mbox{lmgamma}(n, x) =
   \begin{cases}
     \textrm{error} & \mbox{if } x\in \{\dots, -3, -2, -1, 0\}\\
     \ln\Gamma_n(x) & \mbox{if } x\not\in \{\dots, -3, -2, -1, 0\}\\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{lmgamma}(n, x)}{\partial x} =
   \begin{cases}
     \textrm{error} & \mbox{if } x\in \{\dots, -3, -2, -1, 0\}\\
     \frac{\partial\, \ln\Gamma_n(x)}{\partial x} & \mbox{if } x\not\in \{\dots,
 -3, -2, -1, 0\}\\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN} \end{cases}
   \f]

   \f[
   \ln\Gamma_n(x) = \pi^{n(n-1)/4} \, \prod_{j=1}^n \Gamma(x + (1 - j)/2)
   \f]

   \f[
   \frac{\partial \, \ln\Gamma_n(x)}{\partial x} = \sum_{j=1}^n \Psi(x + (1 - j)
 / 2) \f]
 *
 * @tparam T type of scalar
 * @param k Number of dimensions.
 * @param x Function argument.
 * @return Natural log of the multivariate gamma function.
 */
template <typename T, require_arithmetic_t<T>* = nullptr>
inline return_type_t<T> lmgamma(int k, T x) {
  return_type_t<T> result = k * (k - 1) * LOG_PI_OVER_FOUR;

  return result + sum(lgamma(x + (1 - Eigen::ArrayXd::LinSpaced(k, 1, k)) / 2));
}

/**
 * Enables the vectorized application of the natural log of the multivariate
 * gamma function, when the first and/or second arguments are containers.
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @param a First input
 * @param b Second input
 * @return Natural log of the multivariate gamma function applied to the two
 * inputs.
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto lmgamma(const T1& a, const T2& b) {
  return apply_scalar_binary(
      a, b, [&](const auto& c, const auto& d) { return lmgamma(c, d); });
}

}  // namespace math
}  // namespace stan
#endif
