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
          require_all_eigen_vector_t<Vec1, Vec2> * = nullptr,
          require_not_var_t<return_type_t<Vec1, Vec2>> * = nullptr>
inline auto dot_product(const Vec1 &v1, const Vec2 &v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);
  return v1.dot(v2);
}

/**
 * Returns the dot product of the specified arrays.
 *
 * @param v1 First array.
 * @param v2 Second array.
 * @param length Length of both arrays.
 */
template <typename Scalar1, typename Scalar2,
          require_all_stan_scalar_t<Scalar1, Scalar2> * = nullptr,
          require_all_not_var_t<Scalar1, Scalar2> * = nullptr>
inline auto dot_product(const Scalar1 *v1, const Scalar2 *v2, size_t length) {
  return_type_t<Scalar1, Scalar2> result = 0;
  for (size_t i = 0; i < length; i++) {
    result += v1[i] * v2[i];
  }
  return result;
}

/**
 * Returns the dot product of the specified arrays.
 *
 * @param v1 First array.
 * @param v2 Second array.
 * @throw std::domain_error if the vectors are not the same size.
 */
template <typename Scalar1, typename Scalar2, typename Alloc1, typename Alloc2,
          require_all_stan_scalar_t<Scalar1, Scalar2> * = nullptr>
inline return_type_t<Scalar1, Scalar2> dot_product(
    const std::vector<Scalar1, Alloc1> &v1,
    const std::vector<Scalar2, Alloc2> &v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);

  if (v1.size() == 0) {
    return 0.0;
  }

  return dot_product(as_column_vector_or_scalar(v1),
                     as_column_vector_or_scalar(v2));
}

}  // namespace math
}  // namespace stan

#endif
