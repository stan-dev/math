#ifndef STAN_MATH_FWD_MAT_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_FWD_MAT_FUN_VALUE_OF_REC_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/fun/value_of_rec.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * Convert a matrix of fvar to a matrix of doubles.
 *
 * T must implement value_of_rec. See
 * test/math/fwd/mat/fun/value_of_rec.cpp for fvar and var usage.
 *
 * @tparam T Type of Eigen matrix
 * @param[in] x Matrix to be converted
 * @return Matrix of values
 **/
template <typename T, require_eigen_fvar<T>...>
inline auto value_of_rec(T&& x) {
  return value_of_rec(std::forward<decltype(x.val())>(x.val())).eval();
}

}  // namespace math
}  // namespace stan

#endif
