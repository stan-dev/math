#ifndef STAN_MATH_REV_FUN_COLUMNS_DOT_PRODUCT_HPP
#define STAN_MATH_REV_FUN_COLUMNS_DOT_PRODUCT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/prim/meta.hpp>

#include <type_traits>

namespace stan {
namespace math {

template <typename T1, int R1, int C1, typename T2, int R2, int C2,
          typename = require_any_var_t<T1, T2>>
inline Eigen::Matrix<return_type_t<T1, T2>, 1, C1> columns_dot_product(
    const Eigen::Matrix<T1, R1, C1>& v1, const Eigen::Matrix<T2, R2, C2>& v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);
  Eigen::Matrix<var, 1, C1> ret(1, v1.cols());
  for (size_type j = 0; j < v1.cols(); ++j) {
    ret(j) = var(new internal::dot_product_vari<T1, T2>(v1.col(j), v2.col(j)));
  }
  return ret;
}

}  // namespace math
}  // namespace stan
#endif
