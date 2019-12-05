#ifndef STAN_MATH_PRIM_MAT_META_LENGTH_HPP
#define STAN_MATH_PRIM_MAT_META_LENGTH_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {

/** \ingroup type_trait
 * Returns the size of the provided Eigen matrix.
 *
 * @param m a const Eigen matrix
 * @tparam T type of matrix.
 * @return the size of the input matrix
 */
template <typename T, typename = require_eigen_t<T>, typename = void>
size_t length(const T& m) {
  return m.size();
}
}  // namespace stan
#endif
