#ifndef STAN_MATH_PRIM_MAT_FUN_COLUMNS_DOT_PRODUCT_HPP
#define STAN_MATH_PRIM_MAT_FUN_COLUMNS_DOT_PRODUCT_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/arr/err/check_matching_sizes.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of the specified vectors.
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <typename T1, typename T2, enable_if_all_eigen<T1, T2>* = nullptr,
          enable_if_any_not_contains_var<T1, T2>* = nullptr>
inline auto columns_dot_product(const T1& v1, const T2& v2) {
  check_matching_dims("columns_dot_product", "v1", v1, "v2", v2);
  return (v1.transpose() * v2).diagonal();
}

}  // namespace math
}  // namespace stan
#endif
