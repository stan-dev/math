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
 * @tparam MatA type of A (real or complex arithmetic scalars)
 * @tparam MatB type of B (real or complex arithmetic scalars)
 * @tparam MatC type of C (real or complex arithmetic scalars)
 * @param[in] A n x n matrix
 * @param[in] B m x m matrix
 * @param[in] C n x m matrix
 * @return complex n x m solution X
 * @throw std::invalid_argument if A or B is not square, or if the
 *   dimensions of C do not match
 */
template <typename MatA, typename MatB, typename MatC,
          require_all_matrix_t<MatA, MatB, MatC>* = nullptr,
          require_all_bt_arithmetic<MatA, MatB, MatC>* = nullptr>
inline auto sylvester_solver(MatA&& A, MatB&& B, MatC&& C) {
  using complex_mat_t
      = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;

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
  Eigen::ComplexSchur<plain_type_t<MatA>> schur_A{eval(std::forward<MatA>(A))};
  auto&& U = schur_A.matrixU();
  auto&& S = schur_A.matrixT();

  Eigen::ComplexSchur<plain_type_t<MatB>> schur_B{eval(std::forward<MatB>(B))};
  auto&& V = schur_B.matrixU();
  auto&& T = schur_B.matrixT();

  // Y holds D = U^H C V initially, overwritten column by column
  auto&& C_ref = to_ref(std::forward<MatC>(C));
  auto Y = (U.adjoint() * C_ref * V).eval();
  plain_type_t<decltype(S)> Sshift = S;
  for (Eigen::Index i = 0; i < m; ++i) {
    Y.col(i).noalias() -= Y.leftCols(i) * T.col(i).head(i);
    Sshift.diagonal() = S.diagonal().array() + T(i, i);
    Sshift.template triangularView<Eigen::Upper>().solveInPlace(Y.col(i));
  }
  return complex_mat_t(U * Y * V.adjoint());
}

}  // namespace math
}  // namespace stan
#endif
