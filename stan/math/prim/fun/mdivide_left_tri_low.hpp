#ifndef STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_TRI_LOW_HPP
#define STAN_MATH_PRIM_FUN_MDIVIDE_LEFT_TRI_LOW_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/mdivide_left_tri.hpp>

namespace stan {
namespace math {

/**
 * Return the result of left dividing the second argument by the
 * first argument.  Calling <code>mdivide_left_tri_low(A,
 * b)</code> with divisor <code>A</code> and dividend
 * <code>b</code> is more arithmetically stable than calling
 * <code>inv(A) * b</code>.
 *
 * @tparam T1 type of the divisor matrix
 * @tparam T2 type of the dividend matrix
 *
 * @param A divisor, an invertible square matrix
 * @param b dividend, a matrix or vector with the same number of
 *   rows as the divisor has columns
 * @return left division of b by A
 * @throws std::invalid_argument if the divisor is not square or
 *   the dividend does not have the same number of rows as the
 *   divisor has columns.
 */
template <typename T1, typename T2, require_all_eigen_t<T1, T2>* = nullptr,
          require_all_not_eigen_vt<is_fvar, T1, T2>* = nullptr>
inline Eigen::Matrix<return_type_t<T1, T2>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left_tri_low(const T1& A, const T2& b) {
  check_square("mdivide_left_tri_low", "A", A);
  check_multiplicable("mdivide_left_tri_low", "A", A, "b", b);
  if (A.rows() == 0) {
    return {0, b.cols()};
  }

  return mdivide_left_tri<Eigen::Lower>(A, b);
}

template <typename T, require_eigen_t<T>* = nullptr,
          require_not_eigen_vt<is_fvar, T>* = nullptr>
inline plain_type_t<T> mdivide_left_tri_low(const T& A) {
  check_square("mdivide_left_tri_low", "A", A);
  if (A.rows() == 0) {
    return {};
  }

  return mdivide_left_tri<Eigen::Lower>(A);
}

}  // namespace math
}  // namespace stan

#endif
