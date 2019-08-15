#ifndef STAN_MATH_REV_MAT_FUN_COLUMNS_DOT_PRODUCT_HPP
#define STAN_MATH_REV_MAT_FUN_COLUMNS_DOT_PRODUCT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/arr/err/check_matching_sizes.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/rev/mat/fun/dot_product.hpp>
#include <stan/math/prim/meta.hpp>

#include <type_traits>

namespace stan {
namespace math {

template <typename T1, typename T2, enable_if_all_eigen<T1, T2>* = nullptr, enable_if_any_contains_var<T1, T2>* = nullptr>
inline auto columns_dot_product(const T1& v1, const T2& v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);
  typename T1::PlainObject ret(1, v1.cols());
  for (int j = 0; j < v1.cols(); ++j) {
    ret(j) = var(new internal::dot_product_vari<typename T1::Scalar, typename T2::Scalar>(v1.col(j), v2.col(j)));
  }
  return ret;
}

}  // namespace math
}  // namespace stan
#endif
