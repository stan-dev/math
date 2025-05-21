#ifndef STAN_MATH_FWD_FUN_MULTIPLY_HPP
#define STAN_MATH_FWD_FUN_MULTIPLY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <stan/math/prim/fun/multiply.hpp>

namespace stan {
namespace math {

template <typename Mat1, typename Mat2,
          require_all_eigen_vt<is_fvar, Mat1, Mat2>* = nullptr,
          require_vt_same<Mat1, Mat2>* = nullptr,
          require_not_eigen_row_and_col_t<Mat1, Mat2>* = nullptr>
inline auto multiply(const Mat1& m1, const Mat2& m2) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
  return (m1 * m2).eval();
}

template <typename Mat1, typename Mat2,
          require_eigen_vt<is_fvar, Mat1>* = nullptr,
          require_eigen_vt<std::is_floating_point, Mat2>* = nullptr,
          require_not_eigen_row_and_col_t<Mat1, Mat2>* = nullptr>
inline auto multiply(const Mat1& m1, const Mat2& m2) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
  Eigen::Matrix<value_type_t<Mat1>, Mat1::RowsAtCompileTime,
                Mat2::ColsAtCompileTime>
      result(m1.rows(), m2.cols());
  for (size_type i = 0; i < m1.rows(); i++) {
    Eigen::Matrix<value_type_t<Mat1>, 1, Mat1::ColsAtCompileTime> crow
        = m1.row(i);
    for (size_type j = 0; j < m2.cols(); j++) {
      result(i, j) = dot_product(crow, m2.col(j));
    }
  }
  return result;
}

template <typename Mat1, typename Mat2,
          require_eigen_vt<std::is_floating_point, Mat1>* = nullptr,
          require_eigen_vt<is_fvar, Mat2>* = nullptr,
          require_not_eigen_row_and_col_t<Mat1, Mat2>* = nullptr>
inline auto multiply(const Mat1& m1, const Mat2& m2) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
  Eigen::Matrix<value_type_t<Mat2>, Mat1::RowsAtCompileTime,
                Mat2::ColsAtCompileTime>
      result(m1.rows(), m2.cols());
  for (size_type i = 0; i < m1.rows(); i++) {
    Eigen::Matrix<double, 1, Mat1::ColsAtCompileTime> crow = m1.row(i);
    for (size_type j = 0; j < m2.cols(); j++) {
      auto ccol = m2.col(j);
      result(i, j) = dot_product(crow, ccol);
    }
  }
  return result;
}

}  // namespace math
}  // namespace stan
#endif
