#ifndef STAN_MATH_FWD_MAT_FUN_VALUE_OF_HPP
#define STAN_MATH_FWD_MAT_FUN_VALUE_OF_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * Convert a matrix of type T to a matrix of partial types.
 *
 * T must implement value_of. See
 * test/math/fwd/mat/fun/value_of.cpp for fvar and var usage.
 *
 * @tparam T Matrix type
 * @param[in] x Matrix to be converted
 * @return Matrix of values
 **/
template <typename T, require_eigen_fvar<T>...>
inline auto value_of(T&& x) {
  return std::forward<T>(x).val().eval();
}

}  // namespace math
}  // namespace stan

#endif
