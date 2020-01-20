#ifndef STAN_MATH_FWD_FUN_MULTIPLY_HPP
#define STAN_MATH_FWD_FUN_MULTIPLY_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/typedefs.hpp>
#include <stan/math/fwd/fun/dot_product.hpp>

namespace stan {
namespace math {

template <
    typename T1, typename T2, typename = require_all_eigen_vt<is_fvar, T1, T2>,
    typename = require_same_t<typename value_type_t<T1>::Scalar,
                              typename value_type_t<T2>::Scalar>,
    typename
    = require_t<bool_constant<T1::ColsAtCompileTime == T2::RowsAtCompileTime>>,
    typename = require_not_t<
        conjunction<is_eigen_row_vector<T1>, is_eigen_col_vector<T2>>>,
    int = 0>
inline auto multiply(const T1& m1, const T2& m2) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
  return m1 * m2;
}

template <typename T1, typename T2, typename = require_all_eigen_t<T1, T2>,
          typename = require_fvar_t<value_type_t<T1>>,
          typename = require_same_vt<double, T2>,
          typename = require_t<
              bool_constant<T1::ColsAtCompileTime == T2::RowsAtCompileTime>>,
          typename = require_not_t<
              conjunction<is_eigen_row_vector<T1>, is_eigen_col_vector<T2>>>,
          int = 0>
inline auto multiply(const T1& m1, const T2& m2) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
  Eigen::Matrix<value_type_t<T1>, T1::RowsAtCompileTime,
                T2::ColsAtCompileTime>
      result(m1.rows(), m2.cols());
  for (size_type i = 0; i < m1.rows(); i++) {
    Eigen::Matrix<value_type_t<T1>, 1, T1::ColsAtCompileTime> crow
        = m1.row(i);
    for (size_type j = 0; j < m2.cols(); j++) {
      auto ccol = m2.col(j);
      result(i, j) = dot_product(crow, ccol);
    }
  }
  return result;
}

template <typename T1, typename T2, typename = require_all_eigen_t<T1, T2>,
          typename = require_same_vt<double, T1>,
          typename = require_fvar_t<value_type_t<T2>>,
          typename = require_t<
              bool_constant<T1::ColsAtCompileTime == T2::RowsAtCompileTime>>,
          typename = require_not_t<
              conjunction<is_eigen_row_vector<T1>, is_eigen_col_vector<T2>>>,
          char = 0>
inline auto multiply(const T1& m1, const T2& m2) {
  check_multiplicable("multiply", "m1", m1, "m2", m2);
  Eigen::Matrix<value_type_t<T2>, T1::RowsAtCompileTime,
                T2::ColsAtCompileTime>
      result(m1.rows(), m2.cols());
  for (size_type i = 0; i < m1.rows(); i++) {
    Eigen::Matrix<double, 1, T1::ColsAtCompileTime> crow = m1.row(i);
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
