#ifndef STAN_MATH_PRIM_PROB_WISHART_CHOLESKY_LPDF_HPP
#define STAN_MATH_PRIM_PROB_WISHART_CHOLESKY_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/dot_self.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <stan/math/prim/fun/lmgamma.hpp>
#include <stan/math/prim/fun/mdivide_left_tri.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {
/** \ingroup multivar_dists
 * Return the natural logarithm of the unnormalized
 * Wishart density of the specified lower-triangular
 * Cholesky factor variate, positive degrees of freedom,
 * and lower-triangular Cholesky factor of the scale matrix.
 *
 * The scale matrix, `L_S`, must be a lower Cholesky factor
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
  using Eigen::Lower;
  using T_L_Y_ref = ref_type_t<T_y>;
  using T_nu_ref = ref_type_t<T_dof>;
  using T_L_S_ref = ref_type_t<T_scale>;
  using T_return = return_type_t<T_y, T_dof, T_scale>;
  static const char* function = "wishart_cholesky_lpdf";
  Eigen::Index k = L_Y.rows();

  check_greater(function, "Degrees of freedom parameter", nu, k - 1);

  check_square(function, "Cholesky random variable", L_Y);
  check_square(function, "Cholesky scale parameter", L_S);
  check_size_match(function, "side length of random variable", L_Y.rows(),
                   "side length of scale parameter", L_S.rows());

  T_L_Y_ref L_Y_ref = L_Y;
  check_cholesky_factor(function, "Cholesky random variable", L_Y_ref);

  T_L_S_ref L_S_ref = L_S;
  check_cholesky_factor(function, "Cholesky scale matrix", L_S_ref);

  T_nu_ref nu_ref = nu;
  T_return lp(0.0);

  if (include_summand<propto, T_dof>::value) {
    lp += k * LOG_TWO * (1 - 0.5 * nu_ref);
    lp += -lmgamma(k, 0.5 * nu_ref);
  }
  if (include_summand<propto, T_dof, T_scale, T_y>::value) {
    auto L_SinvL_Y = mdivide_left_tri<Eigen::Lower>(L_S_ref, L_Y_ref);
    T_return dot_LSinvLY(0.0);
    Eigen::Matrix<T_return, 1, Eigen::Dynamic> linspaced_rv(k);
    T_return nu_minus_1 = nu_ref - 1;

    for (int i = 0; i < k; i++) {
      dot_LSinvLY += dot_self(L_SinvL_Y.row(i).head(i + 1));
      linspaced_rv(i) = nu_minus_1 - i;
    }
    lp += -0.5 * dot_LSinvLY
          + dot_product(linspaced_rv, log(L_Y_ref.diagonal()))
          - nu_ref * sum(log(L_S_ref.diagonal()));
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
