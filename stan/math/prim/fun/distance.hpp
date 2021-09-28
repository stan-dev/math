#ifndef STAN_MATH_PRIM_FUN_DISTANCE_HPP
#define STAN_MATH_PRIM_FUN_DISTANCE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/abs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/squared_distance.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Returns the distance between two scalars.
 *
 * @tparam T1 type of first scalar.
 * @tparam T2 type of second scalar
 * @param x1 First scalar.
 * @param x2 Second scalar.
 * @return Distance between two scalars
 * @throw std::domain_error If the arguments are not finite.
 */
template <typename T1, typename T2,
          require_all_stan_scalar_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> distance(const T1& x1, const T2& x2) {
  check_finite("distance", "x1", x1);
  check_finite("distance", "x2", x2);
  return abs(x1 - x2);
}

/**
 * Returns the distance between the specified vectors.
 *
 * @tparam T1 type of the first vector (must be derived from \c
 * Eigen::MatrixBase and have one compile time dimension equal to 1)
 * @tparam T2 type of the second vector (must be derived from \c
 * Eigen::MatrixBase and have one compile time dimension equal to 1)
 * @param x1 First vector.
 * @param x2 Second vector.
 * @return Distance between the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size.
 */
template <typename T1, typename T2, require_all_vector_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> distance(const T1& x1, const T2& x2) {
  using std::sqrt;
  check_matching_sizes("distance", "x1", x1, "x2", x2);
  return sqrt(squared_distance(x1, x2));
}

}  // namespace math
}  // namespace stan

#endif
