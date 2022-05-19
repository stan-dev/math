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
 * @tparam Scalar1 type of first scalar.
 * @tparam Scalar2 type of second scalar
 * @param x1 First scalar.
 * @param x2 Second scalar.
 * @return Distance between two scalars
 * @throw std::domain_error If the arguments are not finite.
 */
template <typename Scalar1, typename Scalar2,
          require_all_stan_scalar_t<Scalar1, Scalar2>* = nullptr>
inline return_type_t<Scalar1, Scalar2> distance(const Scalar1& x1,
                                                const Scalar2& x2) {
  check_finite("distance", "x1", x1);
  check_finite("distance", "x2", x2);
  return abs(x1 - x2);
}

/**
 * Returns the distance between the specified vectors.
 *
 * @tparam Vec1 type of the first vector (must be derived from \c
 * Eigen::MatrixBase and have one compile time dimension equal to 1)
 * @tparam Vec2 type of the second vector (must be derived from \c
 * Eigen::MatrixBase and have one compile time dimension equal to 1)
 * @param x1 First vector.
 * @param x2 Second vector.
 * @return Distance between the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size.
 */
template <typename Vec1, typename Vec2,
          require_all_vector_t<Vec1, Vec2>* = nullptr>
inline return_type_t<Vec1, Vec2> distance(const Vec1& x1, const Vec2& x2) {
  using std::sqrt;
  check_matching_sizes("distance", "x1", x1, "x2", x2);
  return sqrt(squared_distance(x1, x2));
}

}  // namespace math
}  // namespace stan

#endif
