#ifndef STAN_MATH_REV_FUN_DIMS_HPP
#define STAN_MATH_REV_FUN_DIMS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/dims.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Pushes dimensions of given argument into given result vector.
 *
 * For `var_value` that is the dimensions of inner `vari_value`.
 * @tparam type in `var_value`
 * @param x argument
 * @param result result
 */
template <typename T>
inline void dims(const var_value<T>& x, std::vector<int>& result) {
  dims(*x.vi_, result);
}

/**
 * Pushes dimensions of given argument into given result vector.
 *
 * For `vari_value` containing Eigen type those are the numbers of rows and
 * columns.
 * @param x argument
 * @param result result
 */
template <typename T>
inline void dims(const vari_value<T>& x, std::vector<int>& result) {
  dims(x.val_, result);
}

}  // namespace math
}  // namespace stan

#endif  // DIMS_HPP
