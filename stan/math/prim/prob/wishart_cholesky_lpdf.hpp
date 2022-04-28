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
 * Return the natural logarithm of the unnormalized Wishart density of the specified
 * lower-triangular Cholesky factor variate, positive degrees of freedom, and lower-triangular
 * Cholesky factor of the scale matrix.
 *
 * The scale matrix, L_S, must be a lower Cholesky factor
 * dimension, k, is implicit
 * nu must be greater than k-1
 *
 * The change of variables from the input positive-definite matrix to
 * the Cholesky factor is given in Theorem 2.1.9 in
 * Muirhead, R. J. (2005).
 * Aspects of Multivariate Statistical Theory. Wiley-Interscience.

 * @tparam T_y Cholesky factor matrix
 * @tparam T_dof scalar degrees of freedom
 * @tparam T_scale Cholesky factor matrix
 * @param L_Y lower triangular Cholesky factor of the inverse covariance matrix
 * @param nu scalar degrees of freedom
 * @param L_S lower triangular Cholesky factor of the scale matrix
 * @return natural logarithm of the Wishart density at L_Y given nu and L_S
 * @throw std::domain_error if nu is not greater than k-1
 * @throw std::domain_error if L_S is not a valid Cholesky factor
 */
template <bool propto, typename T_y, typename T_dof, typename T_scale,
          require_stan_scalar_t<T_dof>* = nullptr,
          require_all_matrix_t<T_y, T_scale>* = nullptr>
return_type_t<T_y, T_dof, T_scale> wishart_cholesky_lpdf(const T_y& L_Y,
                                                         const T_dof& nu,
                                                         const T_scale& L_S) {
  using Eigen::Dynamic;
  using Eigen::Lower;
  using Eigen::Matrix;
  using T_L_Y_ref = ref_type_t<T_y>;
  using T_nu_ref = ref_type_t<T_dof>;
  using T_L_S_ref = ref_type_t<T_scale>;
  static const char* function = "wishart_cholesky_lpdf";
  Eigen::Index k = L_Y.rows();
  check_size_match(function, "Rows of random variable", L_Y.rows(),
                   "columns of scale parameter", L_S.rows());

  T_L_Y_ref L_Y_ref = L_Y;
  T_nu_ref nu_ref = nu;
  T_L_S_ref L_S_ref = L_S;

  check_greater(function, "Degrees of freedom parameter", nu_ref, k - 1);
  check_cholesky_factor(function, "random variable", L_Y_ref);
  check_cholesky_factor(function, "scale parameter", L_S_ref);

  return_type_t<T_y, T_dof, T_scale> lp(0.0);

  if (include_summand<propto, T_dof>::value) {
    lp += k * LOG_TWO * (1 - 0.5 * nu_ref);
  }

  if (include_summand<propto, T_dof>::value) {
    lp -= lmgamma(k, 0.5 * nu_ref);
  }

  if (include_summand<propto, T_dof, T_scale, T_y>::value) {
    auto L_SinvL_Y = mdivide_left_tri<Eigen::Lower>(L_S_ref, L_Y_ref);
    for (int i = 0; i < k; i++) {
      lp -= 0.5 * dot_self(L_SinvL_Y.row(i).head(i + 1))
            + nu_ref * log(L_S_ref.coeff(i, i))
            - (nu_ref - i - 1) * log(L_Y_ref.coeff(i, i));
    }
  }

  return lp;
}

template <typename T_y, typename T_dof, typename T_scale>
inline return_type_t<T_y, T_dof, T_scale> wishart_cholesky_lpdf(
    const T_y& LW, const T_dof& nu, const T_scale& L_S) {
  return wishart_cholesky_lpdf<false>(LW, nu, L_S);
}

}  // namespace math
}  // namespace stan
#endif
