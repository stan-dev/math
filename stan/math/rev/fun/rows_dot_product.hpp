#ifndef STAN_MATH_REV_FUN_ROWS_DOT_PRODUCT_HPP
#define STAN_MATH_REV_FUN_ROWS_DOT_PRODUCT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <type_traits>

namespace stan {
namespace math {

template <typename Mat1, typename Mat2,
          require_all_eigen_t<Mat1, Mat2>* = nullptr,
          require_any_eigen_vt<is_var, Mat1, Mat2>* = nullptr>
inline Eigen::Matrix<var, Mat1::RowsAtCompileTime, 1> rows_dot_product(
    const Mat1& v1, const Mat2& v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);
  Eigen::Matrix<var, Mat1::RowsAtCompileTime, 1> ret(v1.rows(), 1);
  for (size_type j = 0; j < v1.rows(); ++j) {
    ret.coeffRef(j) = dot_product(v1.row(j), v2.row(j));
  }
  return ret;
}

}  // namespace math
}  // namespace stan
#endif
