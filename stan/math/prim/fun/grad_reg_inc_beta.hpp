#ifndef STAN_MATH_PRIM_FUN_GRAD_REG_INC_BETA_HPP
#define STAN_MATH_PRIM_FUN_GRAD_REG_INC_BETA_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/inc_beta.hpp>
#include <stan/math/prim/fun/grad_2F1.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Computes the gradients of the regularized incomplete beta
 * function.  Specifically, this function computes gradients of
 * <code>ibeta(a, b, z)</code>, with respect to the arguments
 * <code>a</code> and <code>b</code>.
 *
 * Uses the equivalence to a hypergeometric function. See
 * http://dlmf.nist.gov/8.17#ii
 *
 * @tparam T type of arguments
 * @param[out] g1 partial derivative of <code>ibeta(a, b, z)</code>
 * with respect to <code>a</code>
 * @param[out] g2 partial derivative of <code>ibeta(a, b,
 * z)</code> with respect to <code>b</code>
 * @param[in] a a
 * @param[in] b b
 * @param[in] z z
 * @param[in] digammaA the value of <code>digamma(a)</code>
 * @param[in] digammaB the value of <code>digamma(b)</code>
 * @param[in] digammaSum the value of <code>digamma(a + b)</code>
 * @param[in] betaAB the value of <code>beta(a, b)</code>
 */
template <typename T>
void grad_reg_inc_beta(T& g1, T& g2, const T& a, const T& b, const T& z,
                       const T& digammaA, const T& digammaB,
                       const T& digammaSum, const T& betaAB) {
  using std::exp;

  T c1 = log(z);
  T c2 = log1m(z);
  T c3 = betaAB * inc_beta(a, b, z);
  T C = exp(a * c1 + b * c2) / a;
  T dF1 = 0;
  T dF2 = 0;
  T dF3 = 0;
  T dFz = 0;
  if (value_of_rec(C)) {
    std::forward_as_tuple(dF1, dF2, dF3, dFz)
        = grad_2F1<true>(a + b, 1.0, a + 1, z);
  }

  T dBda = (c1 - 1.0 / a) * c3 + C * (dF1 + dF3);
  T dBdb = c2 * c3 + C * dF1;

  g1 = (dBda - c3 * (digammaA - digammaSum)) / betaAB;
  g2 = (dBdb - c3 * (digammaB - digammaSum)) / betaAB;
}

}  // namespace math
}  // namespace stan
#endif
