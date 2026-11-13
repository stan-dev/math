#ifndef STAN_MATH_PRIM_FUN_SOLVE_LYAPUNOV_HPP
#define STAN_MATH_PRIM_FUN_SOLVE_LYAPUNOV_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <Eigen/Eigenvalues>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Solve the continuous Lyapunov equation A X + X A^H = C for X.
 *
 * Special case of solve_sylvester with B = A^H. The complex Schur
 * decomposition A = U T U^H also gives A^H = U T^H U^H, so only one
 * decomposition is needed. The equation reduces to T Y + Y T^H = U^H C U,
 * solved column by column from last to first, and X = U Y U^H. A unique
 * solution requires lambda_i(A) + conj(lambda_j(A)) != 0 for all
 * eigenvalues of A.
 *
 * @tparam MatA type of A (real or complex arithmetic scalars)
 * @tparam MatC type of C (real or complex arithmetic scalars)
 * @param[in] A n x n matrix
 * @param[in] C n x n matrix
 * @return complex n x n solution X
 * @throw std::invalid_argument if A is not square, or if the dimensions of C
 *   do not match
 */
template <typename MatA, typename MatC,
          require_all_matrix_t<MatA, MatC>* = nullptr,
          require_all_bt_arithmetic<MatA, MatC>* = nullptr>
inline auto solve_lyapunov(MatA&& A, MatC&& C) {
  using complex_mat_t
      = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;

  check_square("solve_lyapunov", "A", A);
  check_size_match("solve_lyapunov", "rows of A", A.rows(), "rows of C",
                   C.rows());
  check_size_match("solve_lyapunov", "rows of A", A.rows(), "cols of C",
                   C.cols());

  if (unlikely(A.rows() == 0)) {
    return complex_mat_t(0, 0);
  }

  const Eigen::Index n = A.rows();

  // ComplexSchur needs a plain matrix (same pattern as complex_schur_decompose)
  Eigen::ComplexSchur<plain_type_t<MatA>> schur_A{eval(std::forward<MatA>(A))};
  auto&& U = schur_A.matrixU();
  auto&& T = schur_A.matrixT();

  // Y holds D = U^H C U initially; T^H is lower triangular, so Y is
  // overwritten column by column from last to first
  auto&& C_ref = to_ref(std::forward<MatC>(C));
  auto Y = (U.adjoint() * C_ref * U).eval();
  plain_type_t<decltype(T)> Tshift = T;
  for (Eigen::Index j = n - 1; j >= 0; --j) {
    const Eigen::Index k = n - 1 - j;  // already-solved columns
    Y.col(j).noalias() -= Y.rightCols(k) * T.row(j).tail(k).adjoint();
    Tshift.diagonal() = T.diagonal().array() + std::conj(T(j, j));
    Tshift.template triangularView<Eigen::Upper>().solveInPlace(Y.col(j));
  }
  return complex_mat_t(U * Y * U.adjoint());
}

}  // namespace math
}  // namespace stan
#endif
