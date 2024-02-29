#ifndef STAN_MATH_PRIM_PROB_INV_WISHART_LPDF_HPP
#define STAN_MATH_PRIM_PROB_INV_WISHART_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log_determinant_ldlt.hpp>
#include <stan/math/prim/fun/mdivide_left_ldlt.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/lmgamma.hpp>
#include <stan/math/prim/fun/trace.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * The log of the Inverse-Wishart density for the given W, degrees
 * of freedom, and scale matrix.
 *
 * The scale matrix, S, must be k x k, symmetric, and semi-positive
 * definite.
 *
 * \f{eqnarray*}{
 W &\sim& \mbox{\sf{Inv-Wishart}}_{\nu} (S) \\
 \log (p (W \, |\, \nu, S) ) &=& \log \left( \left(2^{\nu k/2} \pi^{k (k-1) /4}
 \prod_{i=1}^k{\Gamma (\frac{\nu + 1 - i}{2})} \right)^{-1} \times \left| S
 \right|^{\nu/2} \left| W \right|^{-(\nu + k + 1) / 2}
 \times \exp (-\frac{1}{2} \mbox{tr} (S W^{-1})) \right) \\
 &=& -\frac{\nu k}{2}\log(2) - \frac{k (k-1)}{4} \log(\pi) - \sum_{i=1}^{k}{\log
 (\Gamma (\frac{\nu+1-i}{2}))}
 +\frac{\nu}{2} \log(\det(S)) - \frac{\nu+k+1}{2}\log (\det(W)) - \frac{1}{2}
 \mbox{tr}(S W^{-1}) \f}
 *
 * @tparam T_y type of scalar
 * @tparam T_dof type of degrees of freedom
 * @tparam T_scale type of scale
 * @param W A scalar matrix
 * @param nu Degrees of freedom
 * @param S The scale matrix
 * @return The log of the Inverse-Wishart density at W given nu and S.
 * @throw std::domain_error if nu is not greater than k-1
 * @throw std::domain_error if S is not square, not symmetric, or not
 * semi-positive definite.
 */
template <bool propto, typename T_y, typename T_dof, typename T_scale>
return_type_t<T_y, T_dof, T_scale> inv_wishart_lpdf(const T_y& W,
                                                    const T_dof& nu,
                                                    const T_scale& S) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using T_W_ref = ref_type_t<T_y>;
  using T_nu_ref = ref_type_t<T_dof>;
  using T_S_ref = ref_type_t<T_scale>;
  static constexpr const char* function = "inv_wishart_lpdf";
  check_size_match(function, "Rows of random variable", W.rows(),
                   "columns of scale parameter", S.rows());
  check_square(function, "random variable", W);
  check_square(function, "scale parameter", S);
  Eigen::Index k = S.rows();
  T_nu_ref nu_ref = nu;
  T_S_ref S_ref = S;
  T_W_ref W_ref = W;
  check_greater(function, "Degrees of freedom parameter", nu_ref, k - 1);
  check_symmetric(function, "random variable", W_ref);
  check_symmetric(function, "scale parameter", S_ref);

  auto ldlt_W = make_ldlt_factor(W_ref);
  check_ldlt_factor(function, "LDLT_Factor of random variable", ldlt_W);
  auto ldlt_S = make_ldlt_factor(S_ref);
  check_ldlt_factor(function, "LDLT_Factor of scale parameter", ldlt_S);

  return_type_t<T_y, T_dof, T_scale> lp(0.0);

  if (include_summand<propto, T_dof>::value) {
    lp -= lmgamma(k, 0.5 * nu_ref);
  }
  if (include_summand<propto, T_dof, T_scale>::value) {
    lp += 0.5 * nu_ref * log_determinant_ldlt(ldlt_S);
  }
  if (include_summand<propto, T_y, T_dof, T_scale>::value) {
    lp -= 0.5 * (nu_ref + k + 1.0) * log_determinant_ldlt(ldlt_W);
  }
  if (include_summand<propto, T_y, T_scale>::value) {
    lp -= 0.5 * trace(mdivide_left_ldlt(ldlt_W, S_ref));
  }
  if (include_summand<propto, T_dof, T_scale>::value) {
    lp -= nu_ref * k * HALF_LOG_TWO;
  }
  return lp;
}

template <typename T_y, typename T_dof, typename T_scale>
inline return_type_t<T_y, T_dof, T_scale> inv_wishart_lpdf(const T_y& W,
                                                           const T_dof& nu,
                                                           const T_scale& S) {
  return inv_wishart_lpdf<false>(W, nu, S);
}

}  // namespace math
}  // namespace stan
#endif
