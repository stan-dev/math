#ifndef STAN_MATH_PRIM_PROB_MATRIX_NORMAL_PREC_RNG_HPP
#define STAN_MATH_PRIM_PROB_MATRIX_NORMAL_PREC_RNG_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

namespace stan {
namespace math {

/** \ingroup multivar_dists
 * Sample from the the matrix normal distribution for the given Mu,
 * Sigma and D where Sigma and D are given as precision matrices, not
 * covariance matrices.
 *
 * @param Mu The mean matrix.
 * @param Sigma The mxm inverse covariance matrix (i.e., the precision
 * matrix) of the rows of y.
 * @param D The nxn inverse covariance matrix (i.e., the precision
 * matrix) of the columns of y.
 * @param rng Pseudo-random number generator.
 * @return A sample from the distribution, of type Matrix<double,
 * Dynamic, Dynamic>.
 * @throw std::invalid_argument if Sigma or D are not square.
 * @throw std::domain_error if Sigma or D are not symmetric, not
 * semi-positive definite, or if they contain infinities or NaNs.
 * @tparam RNG Type of pseudo-random number generator.
 */
template <class RNG>
inline Eigen::MatrixXd matrix_normal_prec_rng(const Eigen::MatrixXd &Mu,
                                              const Eigen::MatrixXd &Sigma,
                                              const Eigen::MatrixXd &D,
                                              RNG &rng) {
  using boost::normal_distribution;
  using boost::variate_generator;
  static const char *function = "matrix_normal_prec_rng";
  check_positive(function, "Sigma rows", Sigma.rows());
  check_finite(function, "Sigma", Sigma);
  check_symmetric(function, "Sigma", Sigma);
  check_positive(function, "D rows", D.rows());
  check_finite(function, "D", D);
  check_symmetric(function, "D", D);
  check_size_match(function, "Rows of location parameter", Mu.rows(),
                   "Rows of Sigma", Sigma.rows());
  check_size_match(function, "Columns of location parameter", Mu.cols(),
                   "Rows of D", D.rows());
  check_finite(function, "Location parameter", Mu);

  Eigen::LDLT<Eigen::MatrixXd> Sigma_ldlt(Sigma);
  // Sigma = PS^T LS DS LS^T PS
  // PS a permutation matrix.
  // LS lower triangular with unit diagonal.
  // DS diagonal.
  Eigen::LDLT<Eigen::MatrixXd> D_ldlt(D);
  // D = PD^T LD DD LD^T PD

  check_pos_semidefinite(function, "Sigma", Sigma_ldlt);
  check_pos_semidefinite(function, "D", D_ldlt);

  // If
  // C ~ N[0, I, I]
  // Then
  // A C B ~ N[0, A A^T, B^T B]
  // So to get
  // Y - Mu ~ N[0, Sigma^(-1), D^(-1)]
  // We need to do
  // Y - Mu = Q^T^(-1) C R^(-1)
  // Where Q^T^(-1) and R^(-1) are such that
  // Q^(-1) Q^(-1)^T = Sigma^(-1)
  // R^(-1)^T R^(-1) = D^(-1)
  // We choose:
  // Q^(-1)^T = PS^T LS^T^(-1) sqrt[DS]^(-1)
  // R^(-1) = sqrt[DD]^(-1) LD^(-1) PD
  // And therefore
  // Y - Mu = (PS^T LS^T^(-1) sqrt[DS]^(-1)) C (sqrt[DD]^(-1) LD^(-1) PD)

  int m = Sigma.rows();
  int n = D.rows();

  variate_generator<RNG &, normal_distribution<>> std_normal_rng(
      rng, normal_distribution<>(0, 1));

  // X = sqrt[DS]^(-1) C sqrt[DD]^(-1)
  // X ~ N[0, DS, DD]
  Eigen::MatrixXd X(m, n);
  Eigen::VectorXd row_stddev
      = Sigma_ldlt.vectorD().array().inverse().sqrt().matrix();
  Eigen::VectorXd col_stddev
      = D_ldlt.vectorD().array().inverse().sqrt().matrix();
  for (int col = 0; col < n; ++col) {
    for (int row = 0; row < m; ++row) {
      double stddev = row_stddev(row) * col_stddev(col);
      // C(row, col) = std_normal_rng();
      X(row, col) = stddev * std_normal_rng();
    }
  }

  // Y - Mu = PS^T (LS^T^(-1) X LD^(-1)) PD
  // Y' = LS^T^(-1) X LD^(-1)
  // Y' = LS^T.solve(X) LD^(-1)
  // Y' = (LD^(-1)^T (LS^T.solve(X))^T)^T
  // Y' = (LD^T.solve((LS^T.solve(X))^T))^T
  // Y = Mu + PS^T Y' PD
  Eigen::MatrixXd Y = Mu
                      + (Sigma_ldlt.transpositionsP().transpose()
                         * (D_ldlt.matrixU().solve(
                                (Sigma_ldlt.matrixU().solve(X)).transpose()))
                               .transpose()
                         * D_ldlt.transpositionsP());

  return Y;
}
}  // namespace math
}  // namespace stan
#endif
