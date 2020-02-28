#ifndef STAN_MATH_PRIM_FUN_SQUARED_DISTANCE_HPP
#define STAN_MATH_PRIM_FUN_SQUARED_DISTANCE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/square.hpp>

namespace stan {
namespace math {

/**
 * Returns the squared distance.
 *
 * @param x1 First scalar.
 * @param x2 Second scalar.
 * @return Squared distance between scalars
 * @throw std::domain_error Any scalar is not finite.
 */
template <typename T1, typename T2,
          require_all_stan_scalar_t<T1, T2>* = nullptr,
          require_all_not_var_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> squared_distance(const T1& x1, const T2& x2) {
  check_finite("squared_distance", "x1", x1);
  check_finite("squared_distance", "x2", x2);
  return square(x1 - x2);
}

/**
 * Returns the squared distance between the specified vectors
 * of the same dimensions.
 *
 * @tparam T1 type of the first vector
 * @tparam T2 type of the second vector
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Squared distance between vectors.
 */
template <typename T1, typename T2,
          require_all_eigen_vector_t<T1, T2>* = nullptr,
          require_all_not_eigen_vt<is_var, T1, T2>* = nullptr>
inline return_type_t<T1, T2> squared_distance(const T1& v1, const T2& v2) {
  check_matching_sizes("squared_distance", "v1", v1, "v2", v2);
  return (as_column_vector_or_scalar(v1) - as_column_vector_or_scalar(v2))
      .squaredNorm();
}

}  // namespace math
}  // namespace stan

#endif
