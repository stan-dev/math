#ifndef STAN_MATH_REV_FUN_SYLVESTER_SOLVER_HPP
#define STAN_MATH_REV_FUN_SYLVESTER_SOLVER_HPP

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
 * Solve the Sylvester equation A X + X B = C for X, reverse mode.
 *
 * Forward pass as in prim. With G the adjoint of X (real part: adjoint of
 * Re X, imaginary part: adjoint of Im X) and Z the solution of
 * A^H Z + Z B^H = G, the gradients are
 *   grad_A = -Z X^H,  grad_B = -X^H Z,  grad_C = Z.
 * Complex inputs receive Re and Im of their gradient; real inputs receive
 * only the real part.
 *
 * @tparam MatA type of A (real or complex, var or double scalars)
 * @tparam MatB type of B (real or complex, var or double scalars)
 * @tparam MatC type of C (real or complex, var or double scalars)
 * @param[in] A n x n matrix
 * @param[in] B m x m matrix
 * @param[in] C n x m matrix
 * @return complex n x m solution X with std::complex<var> entries
 * @throw std::invalid_argument if A or B is not square, or if the
 *   dimensions of C do not match
 */
template <typename MatA, typename MatB, typename MatC,
          require_all_matrix_t<MatA, MatB, MatC>* = nullptr,
          require_any_t<is_var<base_type_t<MatA>>, is_var<base_type_t<MatB>>,
                        is_var<base_type_t<MatC>>>* = nullptr>
inline auto sylvester_solver(const MatA& A, const MatB& B, const MatC& C) {
  using cdouble = std::complex<double>;
  using ret_type
      = Eigen::Matrix<std::complex<var>, Eigen::Dynamic, Eigen::Dynamic>;

  check_square("sylvester_solver", "A", A);
  check_square("sylvester_solver", "B", B);
  check_size_match("sylvester_solver", "rows of A", A.rows(), "rows of C",
                   C.rows());
  check_size_match("sylvester_solver", "rows of B", B.rows(), "cols of C",
                   C.cols());
  if (unlikely(A.rows() == 0 || B.rows() == 0)) {
    return ret_type(A.rows(), B.rows());
  }

  auto arena_A = to_arena(A);
  auto arena_B = to_arena(B);
  auto arena_C = to_arena(C);

  const Eigen::Index m = B.rows();

  // ---- forward pass: identical to prim, on values ----
  Eigen::ComplexSchur<plain_type_t<decltype(value_of(arena_A))>> schur_A(
      value_of(arena_A));
  auto arena_U = to_arena(schur_A.matrixU());
  auto arena_S = to_arena(schur_A.matrixT());

  Eigen::ComplexSchur<plain_type_t<decltype(value_of(arena_B))>> schur_B(
      value_of(arena_B));
  auto arena_V = to_arena(schur_B.matrixU());
  auto arena_T = to_arena(schur_B.matrixT());

  // Y holds D = U^H C V initially, overwritten column by column
  auto Y = (arena_U.adjoint() * value_of(arena_C) * arena_V).eval();
  // deep copy: copying an arena_matrix would alias arena_S
  auto Sshift = arena_S.eval();
  for (Eigen::Index i = 0; i < m; ++i) {
    Y.col(i).noalias() -= Y.leftCols(i) * arena_T.col(i).head(i);
    Sshift.diagonal() = arena_S.diagonal().array() + arena_T(i, i);
    Sshift.template triangularView<Eigen::Upper>().solveInPlace(Y.col(i));
  }
  arena_t<Eigen::MatrixXcd> arena_X = arena_U * Y * arena_V.adjoint();

  auto res = to_arena(arena_X.unaryExpr(
      [](const cdouble& x) { return std::complex<var>(x.real(), x.imag()); }));

  // ---- reverse pass ----
  reverse_pass_callback([arena_A, arena_B, arena_C, arena_U, arena_S, arena_V,
                         arena_T, arena_X, res]() mutable {
    const Eigen::Index rows = arena_S.rows();
    const Eigen::Index cols = arena_T.rows();

    // G = adj(Re X) + i * adj(Im X)
    Eigen::MatrixXcd G(rows, cols);
    for (Eigen::Index c = 0; c < cols; ++c) {
      for (Eigen::Index r = 0; r < rows; ++r) {
        G(r, c) = cdouble(res(r, c).real().adj(), res(r, c).imag().adj());
      }
    }

    // W holds E = U^H G V initially; solve S^H W + W T^H = E,
    // overwritten column by column from last to first
    auto W = (arena_U.adjoint() * G * arena_V).eval();
    auto Sshift_r = arena_S.eval();
    for (Eigen::Index j = cols - 1; j >= 0; --j) {
      const Eigen::Index k = cols - 1 - j;  // already-solved columns
      W.col(j).noalias() -= W.rightCols(k) * arena_T.row(j).tail(k).adjoint();
      Sshift_r.diagonal() = arena_S.diagonal().array() + arena_T(j, j);
      Sshift_r.template triangularView<Eigen::Upper>().adjoint().solveInPlace(
        W.col(j));
    }
    const auto Z = (arena_U * W * arena_V.adjoint()).eval();

    // Re/Im of the gradient go to the adjoints of Re/Im of the input;
    // real inputs only receive Re(grad)
    auto accumulate = [](auto& in, const Eigen::MatrixXcd& g) {
      using in_t = std::decay_t<decltype(in)>;
      if constexpr (is_complex<value_type_t<in_t>>::value) {
        for (Eigen::Index c = 0; c < g.cols(); ++c) {
          for (Eigen::Index r = 0; r < g.rows(); ++r) {
            in(r, c).real().adj() += g(r, c).real();
            in(r, c).imag().adj() += g(r, c).imag();
          }
        }
      } else {
        in.adj() += g.real();
      }
    };

    if constexpr (is_var<base_type_t<MatA>>::value) {
      accumulate(arena_A, -(Z * arena_X.adjoint()));
    }
    if constexpr (is_var<base_type_t<MatB>>::value) {
      accumulate(arena_B, -(arena_X.adjoint() * Z));
    }
    if constexpr (is_var<base_type_t<MatC>>::value) {
      accumulate(arena_C, Z);
    }
  });
  // copy out of the arena: the callback reads adjoints from res, so the
  // caller must not be able to overwrite its entries
  return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif