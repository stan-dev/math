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
template <typename MatA, typename MatB, typename MatC,
          require_all_matrix_t<MatA, MatB, MatC>* = nullptr,
          require_all_t<std::is_arithmetic<base_type_t<MatA>>,
                        std::is_arithmetic<base_type_t<MatB>>,
                        std::is_arithmetic<base_type_t<MatC>>>* = nullptr>
inline auto sylvester_solver(const MatA& A, const MatB& B,
                                         const MatC& C) {
  check_square("sylvester_solver", "A", A);
  check_square("sylvester_solver", "B", B);
  check_size_match("sylvester_solver", "rows of A", A.rows(), "rows of C",
                   C.rows());
  check_size_match("sylvester_solver", "rows of B", B.rows(), "cols of C",
                   C.cols());

  const Eigen::Index m = B.rows();

  Eigen::ComplexSchur<plain_type_t<MatA>> schur_A(A);
  const auto& U = schur_A.matrixU();
  const auto& S = schur_A.matrixT();

  Eigen::ComplexSchur<plain_type_t<MatB>> schur_B(B);
  const auto& V = schur_B.matrixU();
  const auto& T = schur_B.matrixT();

  // Y holds D = U^H C V initially, overwritten column by column
  auto Y = (U.adjoint() * C * V).eval();
  auto Sshift = S;
  for (Eigen::Index i = 0; i < m; ++i) {
    Y.col(i).noalias() -= Y.leftCols(i) * T.col(i).head(i);
    Sshift.diagonal() = S.diagonal().array() + T(i, i);
    Sshift.template triangularView<Eigen::Upper>().solveInPlace(Y.col(i));
  }

  auto X = (U * Y * V.adjoint()).eval();

  return X;
}

}  // namespace math
}  // namespace stan
#endif