#ifndef STAN_MATH_PRIM_FUN_TRANSPOSE_HPP
#define STAN_MATH_PRIM_FUN_TRANSPOSE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Transposes a matrix.
 * @tparam T type of the matrix or expression
 * @param m matrix or expression
 * @return transposed matrix
 */
template <typename T, require_eigen_t<T>* = nullptr>
Eigen::Matrix<value_type_t<T>, T::ColsAtCompileTime,
              T::RowsAtCompileTime> inline transpose(const T& m) {
  return m.transpose();
}

template <typename T, require_stan_scalar_t<T>* = nullptr>
auto&& transpose(T&& x) {
  return std::forward<T>(x);
}
}  // namespace math
}  // namespace stan
#endif
