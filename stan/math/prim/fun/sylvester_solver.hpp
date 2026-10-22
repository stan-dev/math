#ifndef STAN_MATH_PRIM_FUN_SYLVESTER_SOLVER_HPP
#define STAN_MATH_PRIM_FUN_SYLVESTER_SOLVER_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <Eigen/Eigenvalues>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Solve the Sylvester equation A X + X B = C for X.
 *
 * With complex Schur decompositions A = U S U^H and B = V T V^H the
 * equation reduces to S Y + Y T = U^H C V, solved column by column,
 * and X = U Y V^H. A unique solution requires lambda_i(A) + mu_j(B) != 0
 * for all eigenvalues of A and B.
 *
 * @note fvar tangents are lost for exactly-zero subdiagonal entries (e.g.
 *   triangular A or B), as Eigen's Schur deflates on values only. Fix: solve
 *   A Xd + Xd B = Cd - Ad X - X Bd on the values instead of differentiating.
 *
 * @tparam MatA type of A (real or complex scalars)
 * @tparam MatB type of B (real or complex scalars)
 * @tparam MatC type of C (real or complex scalars)
 * @param[in] A n x n matrix
 * @param[in] B m x m matrix
 * @param[in] C n x m matrix
 * @return complex n x m solution X
 * @throw std::invalid_argument if A or B is not square, or if the
 *   dimensions of C do not match
 */
template <
    typename MatA, typename MatB, typename MatC,
    require_all_matrix_t<MatA, MatB, MatC>* = nullptr,
    require_all_not_t<is_var<base_type_t<MatA>>, is_var<base_type_t<MatB>>,
                      is_var<base_type_t<MatC>>>* = nullptr>
inline auto sylvester_solver(const MatA& A, const MatB& B, const MatC& C) {
  // common base scalar: double, or fvar<...> if any input is fvar
  using base_t
      = return_type_t<base_type_t<MatA>, base_type_t<MatB>, base_type_t<MatC>>;
  using complex_t = std::complex<base_t>;
  using complex_mat_t
      = Eigen::Matrix<complex_t, Eigen::Dynamic, Eigen::Dynamic>;
  // A and B are promoted to base_t but stay real or complex as given
  using A_mat_t
      = Eigen::Matrix<std::conditional_t<is_complex<value_type_t<MatA>>::value,
                                         complex_t, base_t>,
                      Eigen::Dynamic, Eigen::Dynamic>;
  using B_mat_t
      = Eigen::Matrix<std::conditional_t<is_complex<value_type_t<MatB>>::value,
                                         complex_t, base_t>,
                      Eigen::Dynamic, Eigen::Dynamic>;

  check_square("sylvester_solver", "A", A);
  check_square("sylvester_solver", "B", B);
  check_size_match("sylvester_solver", "rows of A", A.rows(), "rows of C",
                   C.rows());
  check_size_match("sylvester_solver", "rows of B", B.rows(), "cols of C",
                   C.cols());

  if (unlikely(A.rows() == 0 || B.rows() == 0)) {
    return complex_mat_t(A.rows(), B.rows());
  }

  const Eigen::Index m = B.rows();

  // ComplexSchur needs a plain matrix (same pattern as complex_schur_decompose)
  Eigen::ComplexSchur<A_mat_t> schur_A{
      A_mat_t(A.template cast<typename A_mat_t::Scalar>())};
  const auto& U = schur_A.matrixU();
  const auto& S = schur_A.matrixT();

  Eigen::ComplexSchur<B_mat_t> schur_B{
      B_mat_t(B.template cast<typename B_mat_t::Scalar>())};
  const auto& V = schur_B.matrixU();
  const auto& T = schur_B.matrixT();

  // Y holds D = U^H C V initially, overwritten column by column
  const auto C_c = C.template cast<complex_t>();
  auto Y = (U.adjoint() * C_c * V).eval();
  auto Sshift = S;
  for (Eigen::Index i = 0; i < m; ++i) {
    Y.col(i).noalias() -= Y.leftCols(i) * T.col(i).head(i);
    Sshift.diagonal() = S.diagonal().array() + T(i, i);
    Sshift.template triangularView<Eigen::Upper>().solveInPlace(Y.col(i));
  }

  complex_mat_t X = U * Y * V.adjoint();

  return X;
}

}  // namespace math
}  // namespace stan
#endif
