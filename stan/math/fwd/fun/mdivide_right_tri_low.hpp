#ifndef STAN_MATH_FWD_FUN_MDIVIDE_RIGHT_TRI_LOW_HPP
#define STAN_MATH_FWD_FUN_MDIVIDE_RIGHT_TRI_LOW_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/mdivide_right.hpp>
#include <stan/math/fwd/fun/multiply.hpp>
#include <stan/math/fwd/fun/to_fvar.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/mdivide_right_tri_low.hpp>

namespace stan {
namespace math {

template <typename EigMat1, typename EigMat2,
          require_all_eigen_vt<is_fvar, EigMat1, EigMat2>* = nullptr,
          require_vt_same<EigMat1, EigMat2>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat1>, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_right_tri_low(const EigMat1& A, const EigMat2& b) {
  using T = typename value_type_t<EigMat1>::Scalar;
  constexpr int R1 = EigMat1::RowsAtCompileTime;
  constexpr int C1 = EigMat1::ColsAtCompileTime;
  constexpr int R2 = EigMat2::RowsAtCompileTime;
  constexpr int C2 = EigMat2::ColsAtCompileTime;

  check_square("mdivide_right_tri_low", "b", b);
  check_multiplicable("mdivide_right_tri_low", "A", A, "b", b);
  if (b.size() == 0) {
    return {A.rows(), 0};
  }

  Eigen::Matrix<T, R1, C1> val_A(A.rows(), A.cols());
  Eigen::Matrix<T, R1, C1> deriv_A(A.rows(), A.cols());
  Eigen::Matrix<T, R2, C2> val_b(b.rows(), b.cols());
  Eigen::Matrix<T, R2, C2> deriv_b(b.rows(), b.cols());
  val_b.setZero();
  deriv_b.setZero();

  const Eigen::Ref<const plain_type_t<EigMat1>>& A_ref = A;
  for (int j = 0; j < A.cols(); j++) {
    for (int i = 0; i < A.rows(); i++) {
      val_A.coeffRef(i, j) = A_ref.coeff(i, j).val_;
      deriv_A.coeffRef(i, j) = A_ref.coeff(i, j).d_;
    }
  }

  const Eigen::Ref<const plain_type_t<EigMat2>>& b_ref = b;
  for (int j = 0; j < b.cols(); j++) {
    for (int i = j; i < b.rows(); i++) {
      val_b.coeffRef(i, j) = b_ref.coeff(i, j).val_;
      deriv_b.coeffRef(i, j) = b_ref.coeff(i, j).d_;
    }
  }

  Eigen::Matrix<T, R1, C2> A_mult_inv_b = mdivide_right(val_A, val_b);
  return to_fvar(A_mult_inv_b,
                 mdivide_right(deriv_A, val_b)
                     - A_mult_inv_b * mdivide_right(deriv_b, val_b));
}

template <typename EigMat1, typename EigMat2,
          require_eigen_vt<is_fvar, EigMat1>* = nullptr,
          require_eigen_vt<std::is_arithmetic, EigMat2>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat1>, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_right_tri_low(const EigMat1& A, const EigMat2& b) {
  using T = typename value_type_t<EigMat1>::Scalar;
  constexpr int R1 = EigMat1::RowsAtCompileTime;
  constexpr int C1 = EigMat1::ColsAtCompileTime;

  check_square("mdivide_right_tri_low", "b", b);
  check_multiplicable("mdivide_right_tri_low", "A", A, "b", b);
  if (b.size() == 0) {
    return {A.rows(), 0};
  }

  Eigen::Matrix<T, R1, C1> val_A(A.rows(), A.cols());
  Eigen::Matrix<T, R1, C1> deriv_A(A.rows(), A.cols());

  const Eigen::Ref<const plain_type_t<EigMat1>>& A_ref = A;
  for (int j = 0; j < A.cols(); j++) {
    for (int i = 0; i < A.rows(); i++) {
      val_A.coeffRef(i, j) = A_ref.coeff(i, j).val_;
      deriv_A.coeffRef(i, j) = A_ref.coeff(i, j).d_;
    }
  }

  plain_type_t<EigMat2> val_b = b.template triangularView<Eigen::Lower>();

  return to_fvar(mdivide_right(val_A, val_b), mdivide_right(deriv_A, val_b));
}

template <typename EigMat1, typename EigMat2,
          require_eigen_vt<std::is_arithmetic, EigMat1>* = nullptr,
          require_eigen_vt<is_fvar, EigMat2>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat2>, EigMat1::RowsAtCompileTime,
                     EigMat2::ColsAtCompileTime>
mdivide_right_tri_low(const EigMat1& A, const EigMat2& b) {
  using T = typename value_type_t<EigMat2>::Scalar;
  constexpr int R1 = EigMat1::RowsAtCompileTime;
  constexpr int R2 = EigMat2::RowsAtCompileTime;
  constexpr int C2 = EigMat2::ColsAtCompileTime;
  check_square("mdivide_right_tri_low", "b", b);
  check_multiplicable("mdivide_right_tri_low", "A", A, "b", b);
  if (b.size() == 0) {
    return {A.rows(), 0};
  }

  Eigen::Matrix<T, R1, C2> A_mult_inv_b(A.rows(), b.cols());
  Eigen::Matrix<T, R2, C2> val_b(b.rows(), b.cols());
  Eigen::Matrix<T, R2, C2> deriv_b(b.rows(), b.cols());
  val_b.setZero();
  deriv_b.setZero();

  for (int j = 0; j < b.cols(); j++) {
    for (int i = j; i < b.rows(); i++) {
      val_b(i, j) = b(i, j).val_;
      deriv_b(i, j) = b(i, j).d_;
    }
  }

  A_mult_inv_b = mdivide_right(A, val_b);

  return to_fvar(A_mult_inv_b,
                 -multiply(A_mult_inv_b, mdivide_right(deriv_b, val_b)));
}

}  // namespace math
}  // namespace stan
#endif
