#ifndef STAN_MATH_PRIM_MAT_FUN_DISTANCE_HPP
#define STAN_MATH_PRIM_MAT_FUN_DISTANCE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/fun/squared_distance.hpp>
#include <stan/math/prim/mat/err/check_vector.hpp>
#include <stan/math/prim/arr/err/check_matching_sizes.hpp>

namespace stan {
namespace math {

/**
 * Returns the distance between the specified vectors.
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <typename T1, typename T2, enable_if_all_eigen<T1, T2>* = nullptr>
inline auto distance(const T2& v1, const T2& v2) {
  using std::sqrt;
  check_vector("distance", "v1", v1);
  check_vector("distance", "v2", v2);
  check_matching_sizes("distance", "v1", v1, "v2", v2);
  return sqrt(squared_distance(v1, v2));
}

}  // namespace math
}  // namespace stan
#endif
