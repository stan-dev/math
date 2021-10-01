#ifndef STAN_MATH_PRIM_FUN_LOG_SUM_EXP_SIGNED_HPP
#define STAN_MATH_PRIM_FUN_LOG_SUM_EXP_SIGNED_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/to_vector.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <cmath>
#include <vector>

namespace stan {
namespace math {

/**
 * Return the log of the sum of the exponentiated values of the specified
 * matrix of values.  The matrix may be a full matrix, a vector,
 * or a row vector. Additionally, a matching container of 'signs' indicates
 * whether the exponentiated input should be added or substracted
 *
 * The function is defined as follows to prevent overflow in exponential
 * calculations.
 *
 * \f$\log \sum_{n=1}^N \exp(x_n) = \max(x) + \log \sum_{n=1}^N \exp(x_n -
 * \max(x))\f$.
 *
 * @tparam T1 type of input vector or matrix
 * @tparam T2 type of signs vector or matrix
 * @param[in] x container of specified values
 * @param[in] signs container of signs
 * @return The log of the sum of the exponentiated vector values.
 */
template <typename T1, typename T2,
          require_container_st<std::is_arithmetic, T1>* = nullptr,
          require_container_st<std::is_integral, T2>* = nullptr>
inline auto log_sum_exp_signed(const T1& x, const T2& signs) {
  return apply_vector_unary<T1>::reduce(x, [&](const auto& v) {
    if (v.size() == 0) {
      return NEGATIVE_INFTY;
    }
    const auto& v_ref = to_ref(to_vector(v));
    const auto& signs_ref = to_ref(to_vector(signs));
    const double max = v_ref.cwiseProduct(signs_ref).maxCoeff();
    if (!std::isfinite(max)) {
      return max;
    }
    return max
           + std::log((v_ref.array() - max)
                          .exp()
                          .matrix()
                          .cwiseProduct(signs_ref)
                          .sum());
  });
}

}  // namespace math
}  // namespace stan

#endif
