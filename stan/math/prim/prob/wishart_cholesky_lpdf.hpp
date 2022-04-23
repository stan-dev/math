#ifndef STAN_MATH_PRIM_PROB_WISHART_CHOLESKY_LPDF_HPP
#define STAN_MATH_PRIM_PROB_WISHART_CHOLESKY_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/lmgamma.hpp>
#include <stan/math/prim/fun/mdivide_left_tri.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * The log of the Wishart density for the given Cholesky factor, degrees of freedom,
 * and scale Cholesky factor matrix.
 *
 * The scale matrix, LS, must be a lower Cholesky factor.
 * Dimension, k, is implicit.
 * nu must be greater than k-1
 *
 * The change of variables from the input positive-definite matrix to
 * the Cholesky factor is given in Theorem 2.1.9 in 
 * Muirhead, R. J. (2005). 
 * Aspects of Multivariate Statistical Theory. Wiley-Interscience. 

 * @tparam T_y type of matrix
 * @tparam T_dof type of degrees of freedom
 * @tparam T_scale type of scale
 * @param LY A lower triangular Cholesky factor covariance matrix
 * @param nu Degrees of freedom
 * @param LS The Cholesky factor of the scale matrix
 * @return The log of the Wishart density at LY given nu and LS.
 * @throw std::domain_error if nu is not greater than k-1
 * @throw std::domain_error if S is not square, not symmetric, or not
 * semi-positive definite.
 */
template <bool propto, typename T_y, typename T_dof, typename T_scale,
          require_stan_scalar_t<T_dof>* = nullptr,
          require_all_matrix_t<T_y, T_scale>* = nullptr>
return_type_t<T_y, T_dof, T_scale> wishart_cholesky_lpdf(const T_y& LY, const T_dof& nu,
                                                const T_scale& LS) {
  using Eigen::Dynamic;
  using Eigen::Lower;
  using Eigen::Matrix;
  using T_LY_ref = ref_type_t<T_y>;
  using T_nu_ref = ref_type_t<T_dof>;
  using T_LS_ref = ref_type_t<T_scale>;
  static const char* function = "wishart_cholesky_lpdf";
  Eigen::Index k = LY.rows();
  check_size_match(function, "Rows of random variable", LY.rows(),
                   "columns of scale parameter", LS.rows());

  T_LY_ref LY_ref = LY;
  T_nu_ref nu_ref = nu;
  T_LS_ref LS_ref = LS;

  check_greater(function, "Degrees of freedom parameter", nu_ref, k - 1);
  check_cholesky_factor(function, "random variable", LY_ref);
  check_cholesky_factor(function, "scale parameter", LS_ref);

  return_type_t<T_y, T_dof, T_scale> lp(0.0);

  if (include_summand<propto, T_dof>::value) {
     lp += k * LOG_TWO * (1 - 0.5 * nu_ref);
  }

  if (include_summand<propto, T_dof>::value) {
    lp -= lmgamma(k, 0.5 * nu_ref);
  }

  if (include_summand<propto, T_dof, T_scale, T_y>::value) {
    auto LSinvLY = mdivide_left_tri<Eigen::Lower>(LS_ref, LY_ref);
  
    for (int i = 0; i < k; i++) {
       lp -= 0.5 * dot_self(LSinvLY.row(i).head(i + 1)) + nu_ref * log(LS_ref.coeff(i, i)) - (nu_ref - i - 1) * log(LY_ref.coeff(i, i));
    }
  }

  return lp;
}

template <typename T_y, typename T_dof, typename T_scale>
inline return_type_t<T_y, T_dof, T_scale> wishart_cholesky_lpdf(const T_y& LW,
                                                       const T_dof& nu,
                                                       const T_scale& LS) {
  return wishart_cholesky_lpdf<false>(LW, nu, LS);
}

}  // namespace math
}  // namespace stan
#endif
