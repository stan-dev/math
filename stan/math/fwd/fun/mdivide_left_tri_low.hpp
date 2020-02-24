#ifndef STAN_MATH_FWD_FUN_MDIVIDE_LEFT_TRI_LOW_HPP
#define STAN_MATH_FWD_FUN_MDIVIDE_LEFT_TRI_LOW_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/mdivide_left.hpp>
#include <stan/math/fwd/fun/typedefs.hpp>
#include <stan/math/fwd/fun/to_fvar.hpp>
#include <stan/math/fwd/fun/multiply.hpp>
#include <stan/math/fwd/core.hpp>

namespace stan {
namespace math {

template <typename T1, typename T2,
          require_all_eigen_vt<is_fvar, T1, T2>* = nullptr,
          require_same_vt<T1, T2>* = nullptr>
inline Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left_tri_low(const T1& A, const T2& b) {
  using T = typename value_type_t<T1>::Scalar;
  constexpr int R1 = T1::RowsAtCompileTime;
  constexpr int C1 = T1::ColsAtCompileTime;
  constexpr int R2 = T2::RowsAtCompileTime;
  constexpr int C2 = T2::ColsAtCompileTime;

  check_square("mdivide_left_tri_low", "A", A);
  check_multiplicable("mdivide_left_tri_low", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  Eigen::Matrix<T, R1, C2> inv_A_mult_b(A.rows(), b.cols());
  Eigen::Matrix<T, R1, C2> inv_A_mult_deriv_b(A.rows(), b.cols());
  Eigen::Matrix<T, R1, C1> inv_A_mult_deriv_A(A.rows(), A.cols());
  Eigen::Matrix<T, R1, C1> val_A(A.rows(), A.cols());
  Eigen::Matrix<T, R1, C1> deriv_A(A.rows(), A.cols());
  Eigen::Matrix<T, R2, C2> val_b(b.rows(), b.cols());
  Eigen::Matrix<T, R2, C2> deriv_b(b.rows(), b.cols());
  val_A.setZero();
  deriv_A.setZero();

  const Eigen::Ref<const plain_type_t<T1>>& A_ref = A;
  for (size_type j = 0; j < A.cols(); j++) {
    for (size_type i = j; i < A.rows(); i++) {
      val_A(i, j) = A_ref(i, j).val_;
      deriv_A(i, j) = A_ref(i, j).d_;
    }
  }

  const Eigen::Ref<const plain_type_t<T2>>& b_ref = b;
  for (size_type j = 0; j < b.cols(); j++) {
    for (size_type i = 0; i < b.rows(); i++) {
      val_b(i, j) = b_ref(i, j).val_;
      deriv_b(i, j) = b_ref(i, j).d_;
    }
  }

  inv_A_mult_b = mdivide_left(val_A, val_b);
  inv_A_mult_deriv_b = mdivide_left(val_A, deriv_b);
  inv_A_mult_deriv_A = mdivide_left(val_A, deriv_A);

  Eigen::Matrix<T, R1, C2> deriv(A.rows(), b.cols());
  deriv = inv_A_mult_deriv_b - multiply(inv_A_mult_deriv_A, inv_A_mult_b);

  return to_fvar(inv_A_mult_b, deriv);
}

template <typename T1, typename T2, require_eigen_t<T1>* = nullptr,
          require_same_vt<double, T1>* = nullptr,
          require_eigen_vt<is_fvar, T2>* = nullptr>
inline Eigen::Matrix<value_type_t<T2>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left_tri_low(const T1& A, const T2& b) {
  using T = typename value_type_t<T2>::Scalar;
  constexpr int R1 = T1::RowsAtCompileTime;
  constexpr int C1 = T1::ColsAtCompileTime;
  constexpr int R2 = T2::RowsAtCompileTime;
  constexpr int C2 = T2::ColsAtCompileTime;

  check_square("mdivide_left_tri_low", "A", A);
  check_multiplicable("mdivide_left_tri_low", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  Eigen::Matrix<T, R1, C2> inv_A_mult_b(A.rows(), b.cols());
  Eigen::Matrix<T, R1, C2> inv_A_mult_deriv_b(A.rows(), b.cols());
  Eigen::Matrix<T, R2, C2> val_b(b.rows(), b.cols());
  Eigen::Matrix<T, R2, C2> deriv_b(b.rows(), b.cols());
  Eigen::Matrix<double, R1, C1> val_A(A.rows(), A.cols());
  val_A.setZero();

  const Eigen::Ref<const plain_type_t<T1>>& A_ref = A;
  for (size_type j = 0; j < A.cols(); j++) {
    for (size_type i = j; i < A.rows(); i++) {
      val_A(i, j) = A_ref(i, j);
    }
  }

  const Eigen::Ref<const plain_type_t<T2>>& b_ref = b;
  for (size_type j = 0; j < b.cols(); j++) {
    for (size_type i = 0; i < b.rows(); i++) {
      val_b(i, j) = b_ref(i, j).val_;
      deriv_b(i, j) = b_ref(i, j).d_;
    }
  }

  inv_A_mult_b = mdivide_left(val_A, val_b);
  inv_A_mult_deriv_b = mdivide_left(val_A, deriv_b);

  Eigen::Matrix<T, R1, C2> deriv(A.rows(), b.cols());
  deriv = inv_A_mult_deriv_b;

  return to_fvar(inv_A_mult_b, deriv);
}

template <typename T1, typename T2, require_eigen_vt<is_fvar, T1>* = nullptr,
          require_eigen_t<T2>* = nullptr,
          require_same_vt<double, T2>* = nullptr>
inline Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left_tri_low(const T1& A, const T2& b) {
  using T = typename value_type_t<T1>::Scalar;
  constexpr int R1 = T1::RowsAtCompileTime;
  constexpr int C1 = T1::ColsAtCompileTime;
  constexpr int R2 = T2::RowsAtCompileTime;
  constexpr int C2 = T2::ColsAtCompileTime;

  check_square("mdivide_left_tri_low", "A", A);
  check_multiplicable("mdivide_left_tri_low", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  Eigen::Matrix<T, R1, C2> inv_A_mult_b(A.rows(), b.cols());
  Eigen::Matrix<T, R1, C1> inv_A_mult_deriv_A(A.rows(), A.cols());
  Eigen::Matrix<T, R1, C1> val_A(A.rows(), A.cols());
  Eigen::Matrix<T, R1, C1> deriv_A(A.rows(), A.cols());
  val_A.setZero();
  deriv_A.setZero();

  const Eigen::Ref<const plain_type_t<T1>>& A_ref = A;
  for (size_type j = 0; j < A.cols(); j++) {
    for (size_type i = j; i < A.rows(); i++) {
      val_A(i, j) = A_ref(i, j).val_;
      deriv_A(i, j) = A_ref(i, j).d_;
    }
  }

  inv_A_mult_b = mdivide_left(val_A, b);
  inv_A_mult_deriv_A = mdivide_left(val_A, deriv_A);

  Eigen::Matrix<T, R1, C2> deriv(A.rows(), b.cols());
  deriv = -multiply(inv_A_mult_deriv_A, inv_A_mult_b);

  return to_fvar(inv_A_mult_b, deriv);
}

}  // namespace math
}  // namespace stan
#endif
