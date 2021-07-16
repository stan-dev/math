#ifndef STAN_MATH_REV_FUN_APPEND_ROW_HPP
#define STAN_MATH_REV_FUN_APPEND_ROW_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/append_row.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the result of stacking the rows of the first argument
 * matrix on top of the second argument matrix.
 *
 * Given input types result in following outputs:
 * (matrix, matrix) -> matrix,
 * (matrix, row_vector) -> matrix,
 * (row_vector, matrix) -> matrix,
 * (row_vector, row_vector) -> matrix,
 * (vector, vector) -> vector.
 *
 * @tparam T1 A `var_value` with inner matrix type
 * @tparam T1 A `var_value` with inner matrix type
 *
 * @param A First matrix.
 * @param B Second matrix.
 * @return Result of stacking first matrix on top of second.
 */
template <typename T1, typename T2, require_any_var_matrix_t<T1, T2>* = nullptr>
inline auto append_row(const T1& A, const T2& B) {
  check_size_match("append_row", "columns of A", A.cols(), "columns of B",
                   B.cols());
  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    return make_callback_var(
        append_row(value_of(arena_A), value_of(arena_B)),
        [arena_A, arena_B](auto& vi) mutable {
          arena_A.adj() += vi.adj().topRows(arena_A.rows());
          arena_B.adj() += vi.adj().bottomRows(arena_B.rows());
        });
  } else if (!is_constant<T1>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_A = A;
    return make_callback_var(append_row(value_of(arena_A), value_of(B)),
                             [arena_A](auto& vi) mutable {
                               arena_A.adj()
                                   += vi.adj().topRows(arena_A.rows());
                             });
  } else {
    arena_t<promote_scalar_t<var, T2>> arena_B = B;
    return make_callback_var(append_row(value_of(A), value_of(arena_B)),
                             [arena_B](auto& vi) mutable {
                               arena_B.adj()
                                   += vi.adj().bottomRows(arena_B.rows());
                             });
  }
}

/**
 * Return the result of stacking an scalar on top of the
 * a vector, with the result being a vector.
 *
 * This function applies to (scalar, vector) and returns a vector.
 *
 * @tparam Scal type of the scalar
 * @tparam ColVec A `var_value` with inner column vector type.
 *
 * @param A scalar.
 * @param B vector.
 * @return Result of stacking the scalar on top of the vector.
 */
template <typename Scal, typename ColVec,
          require_stan_scalar_t<Scal>* = nullptr,
          require_t<is_eigen_col_vector<ColVec>>* = nullptr>
inline auto append_row(const Scal& A, const var_value<ColVec>& B) {
  if (!is_constant<Scal>::value && !is_constant<ColVec>::value) {
    var arena_A = A;
    arena_t<promote_scalar_t<var, ColVec>> arena_B = B;
    return make_callback_var(append_row(value_of(arena_A), value_of(arena_B)),
                             [arena_A, arena_B](auto& vi) mutable {
                               arena_A.adj() += vi.adj().coeff(0);
                               arena_B.adj() += vi.adj().tail(arena_B.size());
                             });
  } else if (!is_constant<Scal>::value) {
    var arena_A = A;
    return make_callback_var(
        append_row(value_of(arena_A), value_of(B)),
        [arena_A](auto& vi) mutable { arena_A.adj() += vi.adj().coeff(0); });
  } else {
    arena_t<promote_scalar_t<var, ColVec>> arena_B = B;
    return make_callback_var(append_row(value_of(A), value_of(arena_B)),
                             [arena_B](auto& vi) mutable {
                               arena_B.adj() += vi.adj().tail(arena_B.size());
                             });
  }
}

/**
 * Return the result of stacking a vector on top of the
 * an scalar, with the result being a vector.
 *
 * This function applies to (vector, scalar) and returns a vector.
 *
 * @tparam ColVec a `var_value` with inner column vector type.
 * @tparam Scal type of the scalar
 *
 * @param A vector.
 * @param B scalar.
 * @return Result of stacking the vector on top of the scalar.
 */
template <typename ColVec, typename Scal,
          require_t<is_eigen_col_vector<ColVec>>* = nullptr,
          require_stan_scalar_t<Scal>* = nullptr>
inline auto append_row(const var_value<ColVec>& A, const Scal& B) {
  if (!is_constant<ColVec>::value && !is_constant<Scal>::value) {
    arena_t<promote_scalar_t<var, ColVec>> arena_A = A;
    var arena_B = B;
    return make_callback_var(append_row(value_of(arena_A), value_of(arena_B)),
                             [arena_A, arena_B](auto& vi) mutable {
                               arena_A.adj() += vi.adj().head(arena_A.size());
                               arena_B.adj()
                                   += vi.adj().coeff(vi.adj().size() - 1);
                             });
  } else if (!is_constant<ColVec>::value) {
    arena_t<promote_scalar_t<var, ColVec>> arena_A = A;
    return make_callback_var(append_row(value_of(arena_A), value_of(B)),
                             [arena_A](auto& vi) mutable {
                               arena_A.adj() += vi.adj().head(arena_A.size());
                             });
  } else {
    arena_t<promote_scalar_t<var, Scal>> arena_B = B;
    return make_callback_var(append_row(value_of(A), value_of(arena_B)),
                             [arena_B](auto& vi) mutable {
                               arena_B.adj()
                                   += vi.adj().coeff(vi.adj().size() - 1);
                             });
  }
}

}  // namespace math
}  // namespace stan

#endif
