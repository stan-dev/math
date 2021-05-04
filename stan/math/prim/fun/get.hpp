#ifndef STAN_MATH_PRIM_FUN_GET_HPP
#define STAN_MATH_PRIM_FUN_GET_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <vector>

namespace stan {

/** \ingroup type_trait
 * Returns the provided element. Scalar type overload
 * for the function to retrieve n-th element of a vector,
 * \c Eigen \c Matrix or expression
 *
 * @param x input scalar
 * @param n index of the element to return
 * @return input scalar
 */
template <typename T, typename = require_stan_scalar_t<T>>
inline T get(const T& x, size_t n) {
  return x;
}

/** \ingroup type_trait
 * Returns the n-th element of the provided std::vector.
 *
 * @param x input vector
 * @param n index of the element to return
 * @return n-th element of the input vector
 */
template <typename T>
inline T get(const std::vector<T>& x, size_t n) {
  return x[n];
}

/** \ingroup type_trait
 * Returns the n-th element of the provided Eigen matrix.
 *
 * @param m input \c Eigen \c Matrix or expression
 * @param n index of the element to return
 * @return n-th element of the \c Eigen \c Matrix or expression
 */
template <typename T, typename = require_eigen_t<T>>
inline scalar_type_t<T> get(const T& m, size_t n) {
  return m(static_cast<int>(n));
}

}  // namespace stan
#endif
