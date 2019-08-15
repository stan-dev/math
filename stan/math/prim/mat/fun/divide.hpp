#ifndef STAN_MATH_PRIM_MAT_FUN_DIVIDE_HPP
#define STAN_MATH_PRIM_MAT_FUN_DIVIDE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Return specified matrix divided by specified scalar.
 * @tparam R Row type for matrix.
 * @tparam C Column type for matrix.
 * @param m Matrix.
 * @param c Scalar.
 * @return Matrix divided by scalar.
 */
template <typename T1, typename T2, enable_if_eigen<T1>* = nullptr>
inline auto divide(const T1& m, const T2 c) {
  return m / c;
}

}  // namespace math
}  // namespace stan
#endif
