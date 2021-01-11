#ifndef STAN_MATH_PRIM_PROB_MATRIX_NORMAL_PREC_LPDF_HPP
#define STAN_MATH_PRIM_PROB_MATRIX_NORMAL_PREC_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/log_determinant_ldlt.hpp>
#include <stan/math/prim/fun/subtract.hpp>
#include <stan/math/prim/fun/trace_gen_quad_form.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * The log of the matrix normal density for the given y, mu, Sigma and D
 * where Sigma and D are given as precision matrices, not covariance matrices.
 *
 * @tparam T_y type of scalar
 * @tparam T_Mu type of location
 * @tparam T_Sigma type of Sigma
 * @tparam T_D type of D
 *
 * @param y An mxn matrix.
 * @param Mu The mean matrix.
 * @param Sigma The mxm inverse covariance matrix (i.e., the precision matrix)
 * of the rows of y.
 * @param D The nxn inverse covariance matrix (i.e., the precision matrix) of
 * the columns of y.
 * @return The log of the matrix normal density.
 * @throw std::domain_error if Sigma or D are not square, not symmetric,
 * or not semi-positive definite.
 */
template <bool propto, typename T_y, typename T_Mu, typename T_Sigma,
          typename T_D,
          require_all_matrix_t<T_y, T_Mu, T_Sigma, T_D>* = nullptr>
return_type_t<T_y, T_Mu, T_Sigma, T_D> matrix_normal_prec_lpdf(
    const T_y& y, const T_Mu& Mu, const T_Sigma& Sigma, const T_D& D) {
  static const char* function = "matrix_normal_prec_lpdf";
  check_positive(function, "Sigma rows", Sigma.rows());
  check_finite(function, "Sigma", Sigma);
  check_symmetric(function, "Sigma", Sigma);

  auto ldlt_Sigma = make_ldlt_factor(Sigma);
  check_ldlt_factor(function, "LDLT_Factor of Sigma", ldlt_Sigma);
  check_positive(function, "D rows", D.rows());
  check_finite(function, "D", D);
  check_symmetric(function, "D", D);

  auto ldlt_D = make_ldlt_factor(D);
  check_ldlt_factor(function, "LDLT_Factor of D", ldlt_D);
  check_size_match(function, "Rows of random variable", y.rows(),
                   "Rows of location parameter", Mu.rows());
  check_size_match(function, "Columns of random variable", y.cols(),
                   "Columns of location parameter", Mu.cols());
  check_size_match(function, "Rows of random variable", y.rows(),
                   "Rows of Sigma", Sigma.rows());
  check_size_match(function, "Columns of random variable", y.cols(),
                   "Rows of D", D.rows());
  check_finite(function, "Location parameter", Mu);
  check_finite(function, "Random variable", y);

  return_type_t<T_y, T_Mu, T_Sigma, T_D> lp(0.0);

  if (include_summand<propto>::value) {
    lp += NEG_LOG_SQRT_TWO_PI * y.cols() * y.rows();
  }

  if (include_summand<propto, T_Sigma>::value) {
    lp += log_determinant_ldlt(ldlt_Sigma) * (0.5 * y.rows());
  }

  if (include_summand<propto, T_D>::value) {
    lp += log_determinant_ldlt(ldlt_D) * (0.5 * y.cols());
  }

  if (include_summand<propto, T_y, T_Mu, T_Sigma, T_D>::value) {
    lp -= 0.5 * trace_gen_quad_form(D, Sigma, subtract(y, Mu));
  }
  return lp;
}

template <typename T_y, typename T_Mu, typename T_Sigma, typename T_D,
          require_all_matrix_t<T_y, T_Mu, T_Sigma, T_D>* = nullptr>
return_type_t<T_y, T_Mu, T_Sigma, T_D> matrix_normal_prec_lpdf(
    const T_y& y, const T_Mu& Mu, const T_Sigma& Sigma, const T_D& D) {
  return matrix_normal_prec_lpdf<false>(y, Mu, Sigma, D);
}

}  // namespace math
}  // namespace stan
#endif
