#ifndef STAN_MATH_REV_FUN_SOLVE_LYAPUNOV_HPP
#define STAN_MATH_REV_FUN_SOLVE_LYAPUNOV_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <Eigen/Eigenvalues>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Solve the continuous Lyapunov equation A X + X A^H = C for X, reverse mode.
 *
 * Forward pass as in prim, with a single complex Schur decomposition
 * A = U T U^H. With G the adjoint of X (real part: adjoint of Re X,
 * imaginary part: adjoint of Im X) and Z the solution of A^H Z + Z A = G,
 * the gradients are
 *   grad_A = -(Z X^H + Z^H X),  grad_C = Z,
 * where the second term of grad_A comes from A^H in place of B. The adjoint
 * equation reuses the same Schur factors as T^H W + W T = U^H G U.
 * Complex inputs receive Re and Im of their gradient; real inputs receive
 * only the real part.
 *
 * @tparam MatA type of A (real or complex, var or double scalars)
 * @tparam MatC type of C (real or complex, var or double scalars)
 * @param[in] A n x n matrix
 * @param[in] C n x n matrix
 * @return complex n x n solution X with std::complex<var> entries
 * @throw std::invalid_argument if A is not square, or if the dimensions of C
 *   do not match
 */
template <typename MatA, typename MatC,
          require_all_matrix_t<MatA, MatC>* = nullptr,
          require_any_bt_var<MatA, MatC>* = nullptr>
inline auto solve_lyapunov(MatA&& A, MatC&& C) {
  using cdouble = std::complex<double>;
  using ret_type
      = Eigen::Matrix<std::complex<var>, Eigen::Dynamic, Eigen::Dynamic>;

  check_square("solve_lyapunov", "A", A);
  check_size_match("solve_lyapunov", "rows of A", A.rows(), "rows of C",
                   C.rows());
  check_size_match("solve_lyapunov", "rows of A", A.rows(), "cols of C",
                   C.cols());
  if (unlikely(A.rows() == 0)) {
    return ret_type(0, 0);
  }

  auto arena_A = to_arena(std::forward<MatA>(A));
  auto arena_C = to_arena(std::forward<MatC>(C));

  const Eigen::Index n = A.rows();

  // ---- forward pass: identical to prim, on values ----
  Eigen::ComplexSchur<plain_type_t<decltype(value_of(arena_A))>> schur_A(
      value_of(arena_A));
  auto arena_U = to_arena(schur_A.matrixU());
  auto arena_T = to_arena(schur_A.matrixT());

  // Y holds D = U^H C U initially, overwritten column by column from last
  // to first
  auto Y = (arena_U.adjoint() * value_of(arena_C) * arena_U).eval();
  plain_type_t<decltype(arena_T)> Tshift = arena_T;
  for (Eigen::Index j = n - 1; j >= 0; --j) {
    const Eigen::Index k = n - 1 - j;  // already-solved columns
    Y.col(j).noalias() -= Y.rightCols(k) * arena_T.row(j).tail(k).adjoint();
    Tshift.diagonal() = arena_T.diagonal().array() + std::conj(arena_T(j, j));
    Tshift.template triangularView<Eigen::Upper>().solveInPlace(Y.col(j));
  }
  arena_t<Eigen::MatrixXcd> arena_X = arena_U * Y * arena_U.adjoint();

  auto res = to_arena(arena_X.cast<std::complex<var>>());

  // ---- reverse pass ----
  reverse_pass_callback([arena_A, arena_C, arena_U, arena_T, arena_X,
                         res]() mutable {
    const Eigen::Index n = arena_T.rows();

    // G = adj(Re X) + i * adj(Im X)
    Eigen::MatrixXcd G(n, n);
    for (Eigen::Index c = 0; c < n; ++c) {
      for (Eigen::Index r = 0; r < n; ++r) {
        G(r, c) = cdouble(res(r, c).real().adj(), res(r, c).imag().adj());
      }
    }

    // W holds E = U^H G U initially; solve T^H W + W T = E,
    // overwritten column by column from first to last
    auto W = (arena_U.adjoint() * G * arena_U).eval();
    auto Tshift_r = arena_T.eval();
    for (Eigen::Index j = 0; j < n; ++j) {
      W.col(j).noalias() -= W.leftCols(j) * arena_T.col(j).head(j);
      Tshift_r.diagonal()
          = arena_T.diagonal().array() + std::conj(arena_T(j, j));
      Tshift_r.template triangularView<Eigen::Upper>().adjoint().solveInPlace(
          W.col(j));
    }
    const auto Z = (arena_U * W * arena_U.adjoint()).eval();

    // Re/Im of the gradient go to the adjoints of Re/Im of the input;
    // real inputs only receive Re(grad)
    auto accumulate = [](auto& in, const Eigen::MatrixXcd& g) {
      using in_t = std::decay_t<decltype(in)>;
      if constexpr (is_complex<value_type_t<in_t>>::value) {
        const Eigen::Index rows = g.rows();
        const Eigen::Index cols = g.cols();
        for (Eigen::Index c = 0; c < cols; ++c) {
          for (Eigen::Index r = 0; r < rows; ++r) {
            in(r, c).real().adj() += g(r, c).real();
            in(r, c).imag().adj() += g(r, c).imag();
          }
        }
      } else {
        in.adj() += g.real();
      }
    };

    if constexpr (is_var<base_type_t<MatA>>::value) {
      accumulate(arena_A, -(Z * arena_X.adjoint() + Z.adjoint() * arena_X));
    }
    if constexpr (is_var<base_type_t<MatC>>::value) {
      accumulate(arena_C, Z);
    }
  });
  return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif
