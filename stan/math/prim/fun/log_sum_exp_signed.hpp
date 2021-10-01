#ifndef STAN_MATH_PRIM_FUN_LOG_SUM_EXP_SIGNED_HPP
#define STAN_MATH_PRIM_FUN_LOG_SUM_EXP_SIGNED_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {


/**
 * Return the log of the sum of the exponentiated values of the specified
 * matrix of values.  The matrix may be a full matrix, a vector,
 * a row vector, or a container of these.
 *
 * The function is defined as follows to prevent overflow in exponential
 * calculations.
 *
 * \f$\log \sum_{n=1}^N \exp(x_n) = \max(x) + \log \sum_{n=1}^N \exp(x_n -
 * \max(x))\f$.
 *
 * @tparam T type of input vector or matrix
 * @param[in] x matrix of specified values
 * @return The log of the sum of the exponentiated vector values.
 */
template <typename T1, typename T2, require_container_st<std::is_arithmetic, T1>* = nullptr,
          require_container_st<std::is_integral, T2>* = nullptr>
inline auto log_sum_exp_signed(const T1& x, const T2& signs) {
  return apply_vector_unary<T1>::reduce(x, [&](const auto& v) {
    if (v.size() == 0) {
      return NEGATIVE_INFTY;
    }
    const auto& v_ref = to_ref(v);
    const double max = v_ref.cwiseProduct(signs).maxCoeff();
    if (!std::isfinite(max)) {
      return max;
    }
    return max + std::log((v_ref.array() - max).exp().cwiseProduct(signs).sum());
  });
}

}  // namespace math
}  // namespace stan

#endif
