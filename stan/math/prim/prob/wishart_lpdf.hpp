#ifndef STAN_MATH_PRIM_PROB_WISHART_LPDF_HPP
#define STAN_MATH_PRIM_PROB_WISHART_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/lmgamma.hpp>
#include <stan/math/prim/fun/trace.hpp>
#include <stan/math/prim/fun/log_determinant_ldlt.hpp>
#include <stan/math/prim/fun/mdivide_left_ldlt.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * The log of the Wishart density for the given W, degrees of freedom,
 * and scale matrix.
 *
 * The scale matrix, S, must be k x k, symmetric, and semi-positive definite.
 * Dimension, k, is implicit.
 * nu must be greater than k-1
 *
 * \f{eqnarray*}{
 W &\sim& \mbox{\sf{Wishart}}_{\nu} (S) \\
 \log (p (W \, |\, \nu, S) ) &=& \log \left( \left(2^{\nu k/2} \pi^{k (k-1) /4}
 \prod_{i=1}^k{\Gamma (\frac{\nu + 1 - i}{2})} \right)^{-1} \times \left| S
 \right|^{-\nu/2} \left| W \right|^{(\nu - k - 1) / 2}
 \times \exp (-\frac{1}{2} \mbox{tr} (S^{-1} W)) \right) \\
 &=& -\frac{\nu k}{2}\log(2) - \frac{k (k-1)}{4} \log(\pi) - \sum_{i=1}^{k}{\log
 (\Gamma (\frac{\nu+1-i}{2}))}
 -\frac{\nu}{2} \log(\det(S)) + \frac{\nu-k-1}{2}\log (\det(W)) - \frac{1}{2}
 \mbox{tr} (S^{-1}W) \f}
 *
 * @param W A scalar matrix
 * @param nu Degrees of freedom
 * @param S The scale matrix
 * @return The log of the Wishart density at W given nu and S.
 * @throw std::domain_error if nu is not greater than k-1
 * @throw std::domain_error if S is not square, not symmetric, or not
 semi-positive definite.
 * @tparam T_y Type of scalar.
 * @tparam T_dof Type of degrees of freedom.
 * @tparam T_scale Type of scale.
 */
template <bool propto, typename T_y, typename T_dof, typename T_scale>
return_type_t<T_y, T_dof, T_scale> wishart_lpdf(
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& W,
    const T_dof& nu,
    const Eigen::Matrix<T_scale, Eigen::Dynamic, Eigen::Dynamic>& S) {
  static const char* function = "wishart_lpdf";

  using Eigen::Dynamic;
  using Eigen::Lower;
  using Eigen::Matrix;

  typename index_type<Matrix<T_scale, Dynamic, Dynamic> >::type k = W.rows();
  return_type_t<T_y, T_dof, T_scale> lp(0.0);
  check_greater(function, "Degrees of freedom parameter", nu, k - 1);
  check_square(function, "random variable", W);
  check_square(function, "scale parameter", S);
  check_size_match(function, "Rows of random variable", W.rows(),
                   "columns of scale parameter", S.rows());

  LDLT_factor<T_y, Eigen::Dynamic, Eigen::Dynamic> ldlt_W(W);
  check_ldlt_factor(function, "LDLT_Factor of random variable", ldlt_W);

  LDLT_factor<T_scale, Eigen::Dynamic, Eigen::Dynamic> ldlt_S(S);
  check_ldlt_factor(function, "LDLT_Factor of scale parameter", ldlt_S);

  if (include_summand<propto, T_dof>::value) {
    lp -= nu * k * HALF_LOG_TWO;
  }

  if (include_summand<propto, T_dof>::value) {
    lp -= lmgamma(k, 0.5 * nu);
  }

  if (include_summand<propto, T_dof, T_scale>::value) {
    lp -= 0.5 * nu * log_determinant_ldlt(ldlt_S);
  }

  if (include_summand<propto, T_scale, T_y>::value) {
    Matrix<return_type_t<T_y, T_scale>, Dynamic, Dynamic> Sinv_W(
        mdivide_left_ldlt(ldlt_S, static_cast<Matrix<T_y, Dynamic, Dynamic> >(
                                      W.template selfadjointView<Lower>())));
    lp -= 0.5 * trace(Sinv_W);
  }

  if (include_summand<propto, T_y, T_dof>::value && nu != (k + 1)) {
    lp += 0.5 * (nu - k - 1.0) * log_determinant_ldlt(ldlt_W);
  }
  return lp;
}

template <typename T_y, typename T_dof, typename T_scale>
inline return_type_t<T_y, T_dof, T_scale> wishart_lpdf(
    const Eigen::Matrix<T_y, Eigen::Dynamic, Eigen::Dynamic>& W,
    const T_dof& nu,
    const Eigen::Matrix<T_scale, Eigen::Dynamic, Eigen::Dynamic>& S) {
  return wishart_lpdf<false>(W, nu, S);
}

}  // namespace math
}  // namespace stan
#endif
