#ifndef STAN_MATH_FWD_FUN_MDIVIDE_RIGHT_HPP
#define STAN_MATH_FWD_FUN_MDIVIDE_RIGHT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/mdivide_right.hpp>
#include <stan/math/prim/fun/multiply.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/multiply.hpp>
#include <stan/math/fwd/fun/to_fvar.hpp>
#include <stan/math/fwd/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<fvar<T>, R1, C2> mdivide_right(
    const Eigen::Matrix<fvar<T>, R1, C1> &A,
    const Eigen::Matrix<fvar<T>, R2, C2> &b) {
  check_square("mdivide_right", "b", b);
  check_multiplicable("mdivide_right", "A", A, "b", b);

  Eigen::Matrix<T, R1, C2> A_mult_inv_b(A.rows(), b.cols());
  Eigen::Matrix<T, R1, C2> deriv_A_mult_inv_b(A.rows(), b.cols());
  Eigen::Matrix<T, R2, C2> deriv_b_mult_inv_b(b.rows(), b.cols());
  Eigen::Matrix<T, R1, C1> val_A(A.rows(), A.cols());
  Eigen::Matrix<T, R1, C1> deriv_A(A.rows(), A.cols());
  Eigen::Matrix<T, R2, C2> val_b(b.rows(), b.cols());
  Eigen::Matrix<T, R2, C2> deriv_b(b.rows(), b.cols());

  for (int j = 0; j < A.cols(); j++) {
    for (int i = 0; i < A.rows(); i++) {
      val_A(i, j) = A(i, j).val_;
      deriv_A(i, j) = A(i, j).d_;
    }
  }

  for (int j = 0; j < b.cols(); j++) {
    for (int i = 0; i < b.rows(); i++) {
      val_b(i, j) = b(i, j).val_;
      deriv_b(i, j) = b(i, j).d_;
    }
  }

  A_mult_inv_b = mdivide_right(val_A, val_b);
  deriv_A_mult_inv_b = mdivide_right(deriv_A, val_b);
  deriv_b_mult_inv_b = mdivide_right(deriv_b, val_b);

  Eigen::Matrix<T, R1, C2> deriv(A.rows(), b.cols());
  deriv = deriv_A_mult_inv_b - multiply(A_mult_inv_b, deriv_b_mult_inv_b);

  return to_fvar(A_mult_inv_b, deriv);
}

template <typename T, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<fvar<T>, R1, C2> mdivide_right(
    const Eigen::Matrix<fvar<T>, R1, C1> &A,
    const Eigen::Matrix<double, R2, C2> &b) {
  check_square("mdivide_right", "b", b);
  check_multiplicable("mdivide_right", "A", A, "b", b);

  Eigen::Matrix<T, R2, C2> deriv_b_mult_inv_b(b.rows(), b.cols());
  Eigen::Matrix<T, R1, C1> val_A(A.rows(), A.cols());
  Eigen::Matrix<T, R1, C1> deriv_A(A.rows(), A.cols());

  for (int j = 0; j < A.cols(); j++) {
    for (int i = 0; i < A.rows(); i++) {
      val_A(i, j) = A(i, j).val_;
      deriv_A(i, j) = A(i, j).d_;
    }
  }

  return to_fvar(mdivide_right(val_A, b), mdivide_right(deriv_A, b));
}

template <typename T, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<fvar<T>, R1, C2> mdivide_right(
    const Eigen::Matrix<double, R1, C1> &A,
    const Eigen::Matrix<fvar<T>, R2, C2> &b) {
  check_square("mdivide_right", "b", b);
  check_multiplicable("mdivide_right", "A", A, "b", b);
  Eigen::Matrix<T, R1, C2> A_mult_inv_b(A.rows(), b.cols());
  Eigen::Matrix<T, R2, C2> deriv_b_mult_inv_b(b.rows(), b.cols());
  Eigen::Matrix<T, R2, C2> val_b(b.rows(), b.cols());
  Eigen::Matrix<T, R2, C2> deriv_b(b.rows(), b.cols());

  for (int j = 0; j < b.cols(); j++) {
    for (int i = 0; i < b.rows(); i++) {
      val_b(i, j) = b(i, j).val_;
      deriv_b(i, j) = b(i, j).d_;
    }
  }

  A_mult_inv_b = mdivide_right(A, val_b);
  deriv_b_mult_inv_b = mdivide_right(deriv_b, val_b);

  Eigen::Matrix<T, R1, C2> deriv(A.rows(), b.cols());
  deriv = -multiply(A_mult_inv_b, deriv_b_mult_inv_b);

  return to_fvar(A_mult_inv_b, deriv);
}

}  // namespace math
}  // namespace stan
#endif
