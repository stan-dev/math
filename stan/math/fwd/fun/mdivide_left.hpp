#ifndef STAN_MATH_FWD_FUN_MDIVIDE_LEFT_HPP
#define STAN_MATH_FWD_FUN_MDIVIDE_LEFT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/mdivide_left.hpp>
#include <stan/math/prim/fun/multiply.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/multiply.hpp>
#include <stan/math/fwd/fun/to_fvar.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T1, typename T2,
          require_all_eigen_vt<is_fvar, T1, T2>* = nullptr,
          require_vt_same<T1, T2>* = nullptr>
inline Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left(const T1& A, const T2& b) {
  using T = typename value_type_t<T1>::Scalar;
  constexpr int S1 = T1::RowsAtCompileTime;
  constexpr int C2 = T2::ColsAtCompileTime;

  check_square("mdivide_left", "A", A);
  check_multiplicable("mdivide_left", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  Eigen::Matrix<T, S1, C2> inv_A_mult_b(A.rows(), b.cols());
  Eigen::Matrix<T, S1, C2> inv_A_mult_deriv_b(A.rows(), b.cols());
  Eigen::Matrix<T, S1, S1> inv_A_mult_deriv_A(A.rows(), A.cols());
  Eigen::Matrix<T, S1, S1> val_A(A.rows(), A.cols());
  Eigen::Matrix<T, S1, S1> deriv_A(A.rows(), A.cols());

  const Eigen::Ref<const plain_type_t<T2>>& b_ref = b;
  const Eigen::Ref<const plain_type_t<T1>>& A_ref = A;
  for (int j = 0; j < A.cols(); j++) {
    for (int i = 0; i < A.rows(); i++) {
      val_A(i, j) = A_ref(i, j).val_;
      deriv_A(i, j) = A_ref(i, j).d_;
    }
  }

  inv_A_mult_b = mdivide_left(val_A, b_ref.val());
  inv_A_mult_deriv_b = mdivide_left(val_A, b_ref.d());
  inv_A_mult_deriv_A = mdivide_left(val_A, deriv_A);

  Eigen::Matrix<T, S1, C2> deriv(A.rows(), b.cols());
  deriv = inv_A_mult_deriv_b - multiply(inv_A_mult_deriv_A, inv_A_mult_b);

  return to_fvar(inv_A_mult_b, deriv);
}

template <typename T1, typename T2,
          require_eigen_vt<std::is_arithmetic, T1>* = nullptr,
          require_eigen_vt<is_fvar, T2>* = nullptr>
inline Eigen::Matrix<value_type_t<T2>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left(const T1& A, const T2& b) {
  constexpr int S1 = T1::RowsAtCompileTime;
  constexpr int C2 = T2::ColsAtCompileTime;

  check_square("mdivide_left", "A", A);
  check_multiplicable("mdivide_left", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  const Eigen::Ref<const plain_type_t<T2>>& b_ref = b;

  return to_fvar(mdivide_left(A, b_ref.val()), mdivide_left(A, b_ref.d()));
}

template <typename T1, typename T2, require_eigen_vt<is_fvar, T1>* = nullptr,
          require_eigen_vt<std::is_arithmetic, T2>* = nullptr>
inline Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left(const T1& A, const T2& b) {
  using T = typename value_type_t<T1>::Scalar;
  constexpr int S1 = T1::RowsAtCompileTime;
  constexpr int C2 = T2::ColsAtCompileTime;

  check_square("mdivide_left", "A", A);
  check_multiplicable("mdivide_left", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  Eigen::Matrix<T, S1, C2> inv_A_mult_b(A.rows(), b.cols());
  Eigen::Matrix<T, S1, S1> inv_A_mult_deriv_A(A.rows(), A.cols());
  Eigen::Matrix<T, S1, S1> val_A(A.rows(), A.cols());
  Eigen::Matrix<T, S1, S1> deriv_A(A.rows(), A.cols());

  const Eigen::Ref<const plain_type_t<T1>>& A_ref = A;
  for (int j = 0; j < A.cols(); j++) {
    for (int i = 0; i < A.rows(); i++) {
      val_A(i, j) = A_ref(i, j).val_;
      deriv_A(i, j) = A_ref(i, j).d_;
    }
  }

  inv_A_mult_b = mdivide_left(val_A, b);
  inv_A_mult_deriv_A = mdivide_left(val_A, deriv_A);

  Eigen::Matrix<T, S1, C2> deriv(A.rows(), b.cols());
  deriv = -multiply(inv_A_mult_deriv_A, inv_A_mult_b);

  return to_fvar(inv_A_mult_b, deriv);
}

}  // namespace math
}  // namespace stan
#endif
