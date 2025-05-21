#ifndef STAN_MATH_REV_FUN_APPEND_COL_HPP
#define STAN_MATH_REV_FUN_APPEND_COL_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/append_col.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the result of appending the second argument matrix after the
 * first argument matrix, that is, putting them side by side, with
 * the first matrix followed by the second matrix.
 *
 * Given input types result in following outputs:
 * (matrix, matrix) -> matrix,
 * (matrix, vector) -> matrix,
 * (vector, matrix) -> matrix,
 * (vector, vector) -> matrix,
 * (row vector, row vector) -> row_vector.
 *
 * @tparam T1 A `var_value` with inner matrix type.
 * @tparam T1 A `var_value` with inner matrix type.
 *
 * @param A First matrix.
 * @param B Second matrix.
 * @return Result of appending the first matrix followed by the
 * second matrix side by side.
 */
template <typename T1, typename T2, require_any_var_matrix_t<T1, T2>* = nullptr>
inline auto append_col(const T1& A, const T2& B) {
  check_size_match("append_col", "columns of A", A.rows(), "columns of B",
                   B.rows());
  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    return make_callback_var(
        append_col(value_of(arena_A), value_of(arena_B)),
        [arena_A, arena_B](auto& vi) mutable {
          arena_A.adj() += vi.adj().leftCols(arena_A.cols());
          arena_B.adj() += vi.adj().rightCols(arena_B.cols());
        });
  } else if (!is_constant<T1>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    return make_callback_var(append_col(value_of(arena_A), value_of(B)),
                             [arena_A](auto& vi) mutable {
                               arena_A.adj()
                                   += vi.adj().leftCols(arena_A.cols());
                             });
  } else {
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    return make_callback_var(append_col(value_of(A), value_of(arena_B)),
                             [arena_B](auto& vi) mutable {
                               arena_B.adj()
                                   += vi.adj().rightCols(arena_B.cols());
                             });
  }
}

/**
 * Return the result of stacking an scalar on top of the
 * a row vector, with the result being a row vector.
 *
 * This function applies to (scalar, row vector) and returns a
 * row vector.
 *
 * @tparam Scal type of the scalar
 * @tparam RowVec A `var_value` with an inner type of row vector.
 *
 * @param A scalar.
 * @param B row vector.
 * @return Result of stacking the scalar on top of the row vector.
 */
template <typename Scal, typename RowVec,
          require_stan_scalar_t<Scal>* = nullptr,
          require_t<is_eigen_row_vector<RowVec>>* = nullptr>
inline auto append_col(const Scal& A, const var_value<RowVec>& B) {
  if (!is_constant<Scal>::value && !is_constant<RowVec>::value) {
    var arena_A = A;
    arena_t<promote_scalar_t<var, RowVec>> arena_B = B;
    return make_callback_var(append_col(value_of(arena_A), value_of(arena_B)),
                             [arena_A, arena_B](auto& vi) mutable {
                               arena_A.adj() += vi.adj().coeff(0);
                               arena_B.adj() += vi.adj().tail(arena_B.size());
                             });
  } else if (!is_constant<Scal>::value) {
    var arena_A = A;
    return make_callback_var(
        append_col(value_of(arena_A), value_of(B)),
        [arena_A](auto& vi) mutable { arena_A.adj() += vi.adj().coeff(0); });
  } else {
    arena_t<promote_scalar_t<var, RowVec>> arena_B = B;
    return make_callback_var(append_col(value_of(A), value_of(arena_B)),
                             [arena_B](auto& vi) mutable {
                               arena_B.adj() += vi.adj().tail(arena_B.size());
                             });
  }
}

/**
 * Return the result of stacking a row vector on top of the
 * an scalar, with the result being a row vector.
 *
 * This function applies to (row vector, scalar) and returns a
 * row vector.
 *
 * @tparam RowVec A `var_value` with an inner type of row vector.
 * @tparam Scal type of the scalar
 *
 * @param A row vector.
 * @param B scalar.
 * @return Result of stacking the row vector on top of the scalar.
 */
template <typename RowVec, typename Scal,
          require_t<is_eigen_row_vector<RowVec>>* = nullptr,
          require_stan_scalar_t<Scal>* = nullptr>
inline auto append_col(const var_value<RowVec>& A, const Scal& B) {
  if (!is_constant<RowVec>::value && !is_constant<Scal>::value) {
    arena_t<promote_scalar_t<var, RowVec>> arena_A = A;
    var arena_B = B;
    return make_callback_var(append_col(value_of(arena_A), value_of(arena_B)),
                             [arena_A, arena_B](auto& vi) mutable {
                               arena_A.adj() += vi.adj().head(arena_A.size());
                               arena_B.adj()
                                   += vi.adj().coeff(vi.adj().size() - 1);
                             });
  } else if (!is_constant<RowVec>::value) {
    arena_t<promote_scalar_t<var, RowVec>> arena_A = A;
    return make_callback_var(append_col(value_of(arena_A), value_of(B)),
                             [arena_A](auto& vi) mutable {
                               arena_A.adj() += vi.adj().head(arena_A.size());
                             });
  } else {
    var arena_B = B;
    return make_callback_var(append_col(value_of(A), value_of(arena_B)),
                             [arena_B](auto& vi) mutable {
                               arena_B.adj()
                                   += vi.adj().coeff(vi.adj().size() - 1);
                             });
  }
}

}  // namespace math
}  // namespace stan

#endif
