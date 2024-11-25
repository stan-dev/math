#ifndef STAN_MATH_PRIM_FUN_LOG_ADD_EXP_HPP
#define STAN_MATH_PRIM_FUN_LOG_ADD_EXP_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <stan/math/prim/functor/apply_scalar_binary.hpp>
#include <cmath>
#include <vector>
#include <algorithm>

namespace stan {
namespace math {

/**
 * Calculates the elementwise sum of exponentials without overflow.
 *
 * \f$\log (\exp(a) + \exp(b)) = m + \log(\exp(a-m) + \exp(b-m))\f$,
 *
 * where \f$m = max(a, b)\f$.
 *
 * @tparam T1 type of the first variable
 * @tparam T2 type of the second variable
 * @param a the first variable
 * @param b the second variable
 */

template <typename T1, typename T2,
    require_all_not_st_var<T1, T2>* = nullptr,
    require_all_stan_scalar_t<T1, T2>* = nullptr>
inline return_type_t<T1, T2> log_add_exp(const T2& a, const T1& b) {
    if (a == NEGATIVE_INFTY) {
    return b;
    }
    if (b == NEGATIVE_INFTY) {
    return a;
    }
    if (a == INFTY || b == INFTY) {
    return INFTY;
    }

    const double max_val = std::max(a, b);
    return max_val + std::log(std::exp(a - max_val) + std::exp(b - max_val));
}

/**
 * Calculates the element-wise log sum of exponentials for two containers.
 * For vectors a and b, computes log(exp(a[i]) + exp(b[i])) for each element i.
 * If sizes don't match, uses the smaller size.
 *
 * @tparam T1 type of first container
 * @tparam T2 type of second container
 * @param a First input container
 * @param b Second input container
 * @return Container with element-wise log_add_exp results
 */
template <typename T, require_container_st<std::is_arithmetic, T>* = nullptr>
inline auto log_add_exp(const T& a, const T& b) {
  if (a.size() != b.size()) {
    throw std::invalid_argument("Binary function: size of x ("
                                + std::to_string(a.size()) + ") and size of y ("
                                + std::to_string(b.size())
                                + ") must match in size");
  }

  const size_t min_size = std::min(a.size(), b.size());
  using return_t = return_type_t<T>;

  std::vector<return_t> result(min_size);

  for (size_t i = 0; i < min_size; ++i) {
    if (a[i] == NEGATIVE_INFTY) {
      result[i] = b[i];  // log_add_exp(-∞, b) = b
    } else if (b[i] == NEGATIVE_INFTY) {
      result[i] = a[i];  // log_add_exp(a, -∞) = a
    } else if (a[i] == INFTY || b[i] == INFTY) {
      result[i] = INFTY;  // log_add_exp(∞, b) = ∞
    } else {
      // Log-add-exp trick
      const double max_val = std::max(a[i], b[i]);
      result[i]
          = max_val
            + std::log(std::exp(a[i] - max_val) + std::exp(b[i] - max_val));
    }
  }

  return result;
}

/**
 *  Enables the vectorized application of the log_add_exp function,
 * when the first and/or second arguments are containers.
 *
 * @tparam T1
 * @tparam T2
 * @param a
 * @param b
 * @return auto
 */
template <typename T1, typename T2, require_any_container_t<T1, T2>* = nullptr>
inline auto log_add_exp(const T1& a, const T2& b) {
    return apply_scalar_binary(
        a, b, [](const auto& c, const auto& d) { return log_sum_exp(c, d); }
    );
}

}  // namespace math
}  // namespace stan

#endif
