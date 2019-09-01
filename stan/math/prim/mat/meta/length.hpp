#ifndef STAN_MATH_PRIM_MAT_META_LENGTH_HPP
#define STAN_MATH_PRIM_MAT_META_LENGTH_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/scal/meta/require_generics.hpp>

namespace stan {

/**
 * Returns the size of the provided Eigen matrix.
 *
 * @param m a const Eigen matrix
 * @tparam T type of matrix.
 * @tparam R number of rows in the input matrix.
 * @tparam C number of columns in the input matrix.
 * @return the size of the input matrix
 */
template <typename T, require_eigen<T>...>
auto&& length(T&& m) {
  return std::forward<decltype(m.size())>(m.size());
}
}  // namespace stan
#endif
