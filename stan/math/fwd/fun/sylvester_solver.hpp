#ifndef STAN_MATH_FWD_FUN_SYLVESTER_SOLVER_HPP
#define STAN_MATH_FWD_FUN_SYLVESTER_SOLVER_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/value_of.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/sylvester_solver.hpp>
#include <complex>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Solve the Sylvester equation A X + X B = C for X, forward mode.
 *
 * The value X solves the equation on the values of A, B and C. Differentiating
 * the equation gives the tangent as the solution of another Sylvester equation
 * with the same coefficients,
 *   A X' + X' B = C' - A' X - X B'.
 * Both solves run on the inner scalar type, so no Schur decomposition is taken
 * on fvar values.
 *
 * @tparam MatA type of A (real or complex, fvar or arithmetic scalars)
 * @tparam MatB type of B (real or complex, fvar or arithmetic scalars)
 * @tparam MatC type of C (real or complex, fvar or arithmetic scalars)
 * @param[in] A n x n matrix
 * @param[in] B m x m matrix
 * @param[in] C n x m matrix
 * @return complex n x m solution X with std::complex<fvar> entries
 * @throw std::invalid_argument if A or B is not square, or if the
 *   dimensions of C do not match
 */
template <typename MatA, typename MatB, typename MatC,
          require_all_matrix_t<MatA, MatB, MatC>* = nullptr,
          require_any_bt_fvar<MatA, MatB, MatC>* = nullptr>
inline auto sylvester_solver(const MatA& A, const MatB& B, const MatC& C) {
  using fvar_t
      = return_type_t<base_type_t<MatA>, base_type_t<MatB>, base_type_t<MatC>>;
  using ret_type
      = Eigen::Matrix<std::complex<fvar_t>, Eigen::Dynamic, Eigen::Dynamic>;

  check_square("sylvester_solver", "A", A);
  check_square("sylvester_solver", "B", B);
  check_size_match("sylvester_solver", "rows of A", A.rows(), "rows of C",
                   C.rows());
  check_size_match("sylvester_solver", "rows of B", B.rows(), "cols of C",
                   C.cols());
  if (unlikely(A.rows() == 0 || B.rows() == 0)) {
    return ret_type(A.rows(), B.rows());
  }

  using complex_inner_t = std::complex<typename fvar_t::Scalar>;

  // tangent of an fvar or complex<fvar> matrix, staying real or complex
  auto tangent = [](const auto& x) {
    return x
        .unaryExpr([](const auto& v) {
          if constexpr (is_complex<std::decay_t<decltype(v)>>::value) {
            return std::complex<std::decay_t<decltype(v.real().d_)>>(
                v.real().d_, v.imag().d_);
          } else {
            return v.d_;
          }
        })
        .eval();
  };

  const auto& A_ref = to_ref(A);
  const auto& B_ref = to_ref(B);
  const auto& C_ref = to_ref(C);
  const auto& A_val = to_ref(value_of(A_ref));
  const auto& B_val = to_ref(value_of(B_ref));
  using inner_mat_t
      = Eigen::Matrix<complex_inner_t, Eigen::Dynamic, Eigen::Dynamic>;

  auto X = sylvester_solver(A_val, B_val, value_of(C_ref));

  inner_mat_t rhs = inner_mat_t::Zero(X.rows(), X.cols());
  if constexpr (is_fvar<base_type_t<MatC>>::value) {
    rhs += tangent(C_ref);
  }
  if constexpr (is_fvar<base_type_t<MatA>>::value) {
    rhs -= tangent(A_ref) * X;
  }
  if constexpr (is_fvar<base_type_t<MatB>>::value) {
    rhs -= X * tangent(B_ref);
  }
  inner_mat_t X_d = sylvester_solver(A_val, B_val, rhs);
  const Eigen::Index X_rows = X.rows();
  const Eigen::Index X_cols = X.cols();
  ret_type res(X_rows, X_cols);
  for (Eigen::Index j = 0; j < X_cols; ++j) {
    for (Eigen::Index i = 0; i < X_rows; ++i) {
      res(i, j)
          = std::complex<fvar_t>(fvar_t(X(i, j).real(), X_d(i, j).real()),
                                 fvar_t(X(i, j).imag(), X_d(i, j).imag()));
    }
  }
  return res;
}

}  // namespace math
}  // namespace stan
#endif
