#ifndef STAN_MATH_PRIM_FUN_DOT_PRODUCT_HPP
#define STAN_MATH_PRIM_FUN_DOT_PRODUCT_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the dot product of the specified vectors.
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::invalid_argument If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <typename Vec1, typename Vec2,
          typename = require_all_eigen_vector_t<Vec1, Vec2>,
          typename = require_all_not_eigen_vt<is_var, Vec1, Vec2>>
inline return_type_t<Vec1, Vec2> dot_product(const Vec1 &v1, const Vec2 &v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);
  return v1.dot(v2);
}

}  // namespace math
}  // namespace stan

#endif
