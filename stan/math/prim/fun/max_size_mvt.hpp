#ifndef STAN_MATH_PRIM_FUN_MAX_SIZE_MVT_HPP
#define STAN_MATH_PRIM_FUN_MAX_SIZE_MVT_HPP

#include <stan/math/prim/fun/size_mvt.hpp>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Calculate the size of the largest multivariate input. A multivariate
 * container is either an Eigen matrix, whose mvt size is 1, or an std::vector
 * of Eigen matrices, whose mvt size is the size of the std::vector. It is an
 * error to supply any other type of input.
 * @tparam T1 type of the first input
 * @tparam Ts types of the other inputs
 * @param x1 first input
 * @param xs other inputs
 * @return the size of the largest input
 * @throw `invalid_argument` if provided with an input that is not an Eigen
 * matrix or std::vector of Eigen matrices
 */
template <typename T1, typename... Ts>
inline size_t max_size_mvt(const T1& x1, const Ts&... xs) {
  return std::max({stan::math::size_mvt(x1), stan::math::size_mvt(xs)...});
}

}  // namespace math
}  // namespace stan
#endif
