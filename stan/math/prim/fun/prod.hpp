#ifndef STAN_MATH_PRIM_FUN_PROD_HPP
#define STAN_MATH_PRIM_FUN_PROD_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the product of given scalar. This is a no-op.
 *
 * @tparam T type of the scalar
 * @param v specified scalar
 * @return the scalar
 */
template <typename T, require_stan_scalar_t<T>* = nullptr>
T prod(const T& v) {
  return v;
}

/**
 * Returns the product of the coefficients of the specified
 * standard vector.
 *
 * @tparam T type of elements in the vector
 * @param v Specified vector.
 * @return Product of coefficients of vector.
 */
template <typename T>
inline T prod(const std::vector<T>& v) {
  if (v.size() == 0) {
    return 1;
  }
  Eigen::Map<const Eigen::Matrix<T, Eigen::Dynamic, 1>> m(&v[0], v.size());
  return m.prod();
}

/**
 * Returns the product of the coefficients of the specified
 * column vector.
 *
 * @tparam T type of elements in the vector
 * @param v Specified vector.
 * @return Product of coefficients of vector.
 */
template <typename EigMat, require_eigen_t<EigMat>* = nullptr>
inline value_type_t<EigMat> prod(const EigMat& v) {
  if (v.size() == 0) {
    return 1.0;
  }
  return v.prod();
}

}  // namespace math
}  // namespace stan

#endif
