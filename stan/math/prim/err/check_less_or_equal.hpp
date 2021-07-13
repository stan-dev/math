#ifndef STAN_MATH_PRIM_ERR_CHECK_LESS_OR_EQUAL_HPP
#define STAN_MATH_PRIM_ERR_CHECK_LESS_OR_EQUAL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/throw_domain_error_vec.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <string>

namespace stan {
namespace math {

/**
 * Check if <code>y</code> is less or equal to <code>high</code>.
 * This function is vectorized and will check each element of
 * <code>y</code> against each element of <code>high</code>.
 * @tparam T_y Type of y
 * @tparam T_high Type of upper bound
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw <code>std::domain_error</code> if y is not less than or equal
 *   to low or if any element of y or high is NaN
 */
template <typename T_y, typename T_high, require_container_t<T_y>* = nullptr,
          require_not_var_matrix_t<T_high>* = nullptr>
inline void check_less_or_equal(const char* function, const char* name,
                                const T_y& y, const T_high& high) {
  const auto& high_ref = to_ref(value_of(high));
  scalar_seq_view<decltype(high_ref)> high_vec(high_ref);
  const auto& y_ref = to_ref(value_of(y));
  for (size_t n = 0; n < stan::math::size(y_ref); n++) {
    if (!(stan::get(y_ref, n) <= high_vec[n])) {
      [&]() STAN_COLD_PATH {
        std::stringstream msg;
        msg << ", but must be less than or equal to ";
        msg << high_vec[n];
        std::string msg_str(msg.str());
        throw_domain_error_vec(function, name, y_ref, n, "is ",
                               msg_str.c_str());
      }();
    }
  }
}

/**
 * Check if <code>y</code> is less or equal to <code>high</code>.
 * This function is vectorized and will check each element of
 * <code>y</code> against each element of <code>high</code>.
 * @tparam T_y Type of y
 * @tparam T_high Type of upper bound
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param high Upper bound
 * @throw <code>std::domain_error</code> if y is not less than or equal
 *   to low or if any element of y or high is NaN
 */
template <typename T_y, typename T_high, require_stan_scalar_t<T_y>* = nullptr,
          require_not_var_matrix_t<T_high>* = nullptr>
inline void check_less_or_equal(const char* function, const char* name,
                                const T_y& y, const T_high& high) {
  const auto& high_ref = to_ref(value_of(high));
  scalar_seq_view<decltype(high_ref)> high_vec(high_ref);
  for (size_t n = 0; n < stan::math::size(high); n++) {
    if (!(y <= high_vec[n])) {
      [&]() STAN_COLD_PATH {
        std::stringstream msg;
        msg << ", but must be less than or equal to ";
        msg << high_vec[n];
        std::string msg_str(msg.str());
        throw_domain_error(function, name, y, "is ", msg_str.c_str());
      }();
    }
  }
}

template <typename T_y, typename T_high,
          require_any_var_matrix_t<T_y, T_high>* = nullptr>
inline void check_less_or_equal(const char* function, const char* name,
                                const T_y& y, const T_high& high) {
  check_less_or_equal(function, name, value_of(y), value_of(high));
}
}  // namespace math
}  // namespace stan
#endif
