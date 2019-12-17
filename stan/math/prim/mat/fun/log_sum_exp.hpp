#ifndef STAN_MATH_PRIM_MAT_FUN_LOG_SUM_EXP_HPP
#define STAN_MATH_PRIM_MAT_FUN_LOG_SUM_EXP_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/vectorize/apply_vector_unary.hpp>
#include <vector>
#include <cmath>
#include <limits>

namespace stan {
namespace math {

/**
 * Return the log of the sum of the exponentiated values of the specified
 * matrix of values.  The matrix may be a full matrix, a vector,
 * or a row vector.
 *
 * The function is defined as follows to prevent overflow in exponential
 * calculations.
 *
 * \f$\log \sum_{n=1}^N \exp(x_n) = \max(x) + \log \sum_{n=1}^N \exp(x_n -
 * \max(x))\f$.
 *
 * @param[in] x Matrix of specified values
 * @return The log of the sum of the exponentiated vector values.
 */
template <typename T, require_t<std::is_arithmetic<scalar_type_t<T>>>...>
inline auto log_sum_exp(T&& x) {
  return apply_vector_unary<T>::reduce(std::forward<T>(x), [](auto& v){
    if (v.size() == 0) {
      return -std::numeric_limits<double>::infinity();
    }

    const double max = v.maxCoeff();
    if (!std::isfinite(max)) {
      return max;
    }
    return max + std::log((v.array() - max).exp().sum());
  });
}

}  // namespace math
}  // namespace stan

#endif
