#ifndef STAN_MATH_FWD_FUN_COLUMNS_DOT_PRODUCT_HPP
#define STAN_MATH_FWD_FUN_COLUMNS_DOT_PRODUCT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/fun/typedefs.hpp>
#include <stan/math/fwd/core.hpp>

namespace stan {
namespace math {

template <typename T, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<fvar<T>, 1, C1> columns_dot_product(
    const Eigen::Matrix<fvar<T>, R1, C1>& v1,
    const Eigen::Matrix<fvar<T>, R2, C2>& v2) {
  check_matching_dims("columns_dot_product", "v1", v1, "v2", v2);
  Eigen::Matrix<fvar<T>, 1, C1> ret(1, v1.cols());
  for (size_type j = 0; j < v1.cols(); ++j) {
    Eigen::Matrix<fvar<T>, R1, C1> ccol1 = v1.col(j);
    Eigen::Matrix<fvar<T>, R2, C2> ccol2 = v2.col(j);
    ret(0, j) = dot_product(ccol1, ccol2);
  }
  return ret;
}

template <typename T, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<fvar<T>, 1, C1> columns_dot_product(
    const Eigen::Matrix<fvar<T>, R1, C1>& v1,
    const Eigen::Matrix<double, R2, C2>& v2) {
  check_matching_dims("columns_dot_product", "v1", v1, "v2", v2);
  Eigen::Matrix<fvar<T>, 1, C1> ret(1, v1.cols());
  for (size_type j = 0; j < v1.cols(); ++j) {
    Eigen::Matrix<fvar<T>, R1, C1> ccol1 = v1.col(j);
    Eigen::Matrix<double, R2, C2> ccol = v2.col(j);
    ret(0, j) = dot_product(ccol1, ccol);
  }
  return ret;
}

template <typename T, int R1, int C1, int R2, int C2>
inline Eigen::Matrix<fvar<T>, 1, C1> columns_dot_product(
    const Eigen::Matrix<double, R1, C1>& v1,
    const Eigen::Matrix<fvar<T>, R2, C2>& v2) {
  check_matching_dims("columns_dot_product", "v1", v1, "v2", v2);
  Eigen::Matrix<fvar<T>, 1, C1> ret(1, v1.cols());
  for (size_type j = 0; j < v1.cols(); ++j) {
    Eigen::Matrix<double, R1, C1> ccol = v1.col(j);
    Eigen::Matrix<fvar<T>, R2, C2> ccol2 = v2.col(j);
    ret(0, j) = dot_product(ccol, ccol2);
  }
  return ret;
}

}  // namespace math
}  // namespace stan
#endif
