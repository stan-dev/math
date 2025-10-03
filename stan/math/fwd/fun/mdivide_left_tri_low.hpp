#ifndef STAN_MATH_FWD_FUN_MDIVIDE_LEFT_TRI_LOW_HPP
#define STAN_MATH_FWD_FUN_MDIVIDE_LEFT_TRI_LOW_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/mdivide_left.hpp>
#include <stan/math/fwd/fun/to_fvar.hpp>
#include <stan/math/fwd/fun/typedefs.hpp>
#include <stan/math/fwd/fun/multiply.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/mdivide_left_tri_low.hpp>

namespace stan {
namespace math {

template <typename T1, typename T2,
          require_all_eigen_vt<is_fvar, T1, T2>* = nullptr,
          require_vt_same<T1, T2>* = nullptr>
inline Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left_tri_low(const T1& A, const T2& b) {
  using T = typename value_type_t<T1>::Scalar;
  constexpr int S1 = T1::RowsAtCompileTime;
  constexpr int C2 = T2::ColsAtCompileTime;

  check_square("mdivide_left_tri_low", "A", A);
  check_multiplicable("mdivide_left_tri_low", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  Eigen::Matrix<T, S1, S1> val_A(A.rows(), A.cols());
  Eigen::Matrix<T, S1, S1> deriv_A(A.rows(), A.cols());
  val_A.setZero();
  deriv_A.setZero();

  const Eigen::Ref<const plain_type_t<T2>>& b_ref = b;
  const Eigen::Ref<const plain_type_t<T1>>& A_ref = A;
  for (size_type j = 0; j < A.cols(); j++) {
    for (size_type i = j; i < A.rows(); i++) {
      val_A(i, j) = A_ref(i, j).val_;
      deriv_A(i, j) = A_ref(i, j).d_;
    }
  }

  Eigen::Matrix<T, S1, C2> inv_A_mult_b = mdivide_left(val_A, b_ref.val());

  return to_fvar(inv_A_mult_b,
                 mdivide_left(val_A, b_ref.d())
                     - multiply(mdivide_left(val_A, deriv_A), inv_A_mult_b));
}

template <typename T1, typename T2, require_eigen_t<T1>* = nullptr,
          require_vt_same<double, T1>* = nullptr,
          require_eigen_vt<is_fvar, T2>* = nullptr>
inline Eigen::Matrix<value_type_t<T2>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left_tri_low(const T1& A, const T2& b) {
  constexpr int S1 = T1::RowsAtCompileTime;

  check_square("mdivide_left_tri_low", "A", A);
  check_multiplicable("mdivide_left_tri_low", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  Eigen::Matrix<double, S1, S1> val_A(A.rows(), A.cols());
  val_A.setZero();

  const Eigen::Ref<const plain_type_t<T2>>& b_ref = b;
  const Eigen::Ref<const plain_type_t<T1>>& A_ref = A;
  for (size_type j = 0; j < A.cols(); j++) {
    for (size_type i = j; i < A.rows(); i++) {
      val_A(i, j) = A_ref(i, j);
    }
  }

  return to_fvar(mdivide_left(val_A, b_ref.val()),
                 mdivide_left(val_A, b_ref.d()));
}

template <typename T1, typename T2, require_eigen_vt<is_fvar, T1>* = nullptr,
          require_eigen_t<T2>* = nullptr,
          require_vt_same<double, T2>* = nullptr>
inline Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime,
                     T2::ColsAtCompileTime>
mdivide_left_tri_low(const T1& A, const T2& b) {
  using T = typename value_type_t<T1>::Scalar;
  constexpr int S1 = T1::RowsAtCompileTime;
  constexpr int C2 = T2::ColsAtCompileTime;

  check_square("mdivide_left_tri_low", "A", A);
  check_multiplicable("mdivide_left_tri_low", "A", A, "b", b);
  if (A.size() == 0) {
    return {0, b.cols()};
  }

  Eigen::Matrix<T, S1, S1> val_A(A.rows(), A.cols());
  Eigen::Matrix<T, S1, S1> deriv_A(A.rows(), A.cols());
  val_A.setZero();
  deriv_A.setZero();

  const Eigen::Ref<const plain_type_t<T1>>& A_ref = A;
  for (size_type j = 0; j < A.cols(); j++) {
    for (size_type i = j; i < A.rows(); i++) {
      val_A(i, j) = A_ref(i, j).val_;
      deriv_A(i, j) = A_ref(i, j).d_;
    }
  }

  Eigen::Matrix<T, S1, C2> inv_A_mult_b = mdivide_left(val_A, b);

  return to_fvar(inv_A_mult_b,
                 -multiply(mdivide_left(val_A, deriv_A), inv_A_mult_b));
}

}  // namespace math
}  // namespace stan
#endif
