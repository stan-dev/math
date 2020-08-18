#ifndef STAN_MATH_PRIM_FUN_MEAN_HPP
#define STAN_MATH_PRIM_FUN_MEAN_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/functor/apply_vector_unary.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Returns the sample mean (i.e., average) of the coefficients
 * in the specified std vector, vector, row vector, or matrix.
 *
 * @tparam T type of the matrix
 *
 * @param m Specified std vector, vector, row vector, or matrix.
 * @return Sample mean of container coefficients.
 */
template <typename T, require_container_t<T>* = nullptr>
inline return_type_t<T> mean(const T& m) {
  check_nonzero_size("mean", "m", m);
  return apply_vector_unary<T>::reduce(m,
                                       [](const auto& a) { return a.mean(); });
}

}  // namespace math
}  // namespace stan

#endif
