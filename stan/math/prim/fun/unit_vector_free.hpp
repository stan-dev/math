#ifndef STAN_MATH_PRIM_FUN_UNIT_VECTOR_FREE_HPP
#define STAN_MATH_PRIM_FUN_UNIT_VECTOR_FREE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Transformation of a unit length vector to a "free" vector
 * However, we are just fixing the unidentified radius to 1.
 * Thus, the transformation is just the identity
 *
 * @tparam T type of elements in the vector
 * @param x unit vector of dimension K
 * @return Unit vector of dimension K considered "free"
 */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1> unit_vector_free(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& x) {
  check_unit_vector("stan::math::unit_vector_free", "Unit vector variable", x);
  return x;
}

}  // namespace math
}  // namespace stan

#endif
