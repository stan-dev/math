#ifndef STAN_MATH_PRIM_MAT_FUN_DOT_PRODUCT_HPP
#define STAN_MATH_PRIM_MAT_FUN_DOT_PRODUCT_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/err/check_vector.hpp>
#include <stan/math/prim/arr/err/check_matching_sizes.hpp>
#include <vector>

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
template <typename T1, typename T2, enable_if_any_eigen<T1, T2> * = nullptr>
inline auto dot_product(const T1 &v1, const T2 &v2) {
  check_vector("dot_product", "v1", v1);
  check_vector("dot_product", "v2", v2);
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);
  return v1.dot(v2);
}

/**
 * Returns the dot product of the specified vectors.
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @param length Number of elements of vectors to sum.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <typename T1, typename T2, typename T3,
          enable_if_all_eigen<T1, T2> * = nullptr>
inline auto dot_product(const T1 &v1, const T2 &v2, const T3 length) {
  return v1.head(length).dot(v2.head(length));
}

/**
 * Returns the dot product of the specified arrays of doubles.
 * @param v1 First array.
 * @param v2 Second array.
 * @param length Number of elements of vectors to sum.
 * @param length Length of both arrays.
 */
template <typename T1, typename T2,
          enable_if_all_contains_stan_scalar<T1, T2> * = nullptr>
inline auto dot_product(const T1 *v1, const T2 *v2, size_t length) {
  return_type_t<T1, T2> result = 0;
  for (size_t i = 0; i < length; i++)
    result += v1[i] * v2[i];
  return result;
}
/**
 * Returns the dot product of the specified arrays of doubles.
 * @param v1 First array.
 * @param v2 Second array.
 * @throw std::domain_error if the vectors are not the same size.
 */
template <typename T1, typename T2>
inline auto dot_product(const std::vector<T1> &v1, const std::vector<T2> &v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);
  return dot_product(&v1[0], &v2[0], v1.size());
}

// Brought in from fvar imp
// This is super weird tho???
template <typename T1, typename T2, typename T3>
inline auto dot_product(const std::vector<T1> &v1, const std::vector<T2> &v2,
                        const T3 length) {
  return dot_product(&v1[0], &v2[0], length);
}
}  // namespace math
}  // namespace stan
#endif
