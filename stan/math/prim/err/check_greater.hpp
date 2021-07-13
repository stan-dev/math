#ifndef STAN_MATH_PRIM_ERR_CHECK_GREATER_HPP
#define STAN_MATH_PRIM_ERR_CHECK_GREATER_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/err/throw_domain_error_vec.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <string>

namespace stan {
namespace math {

/**
 * Check if <code>y</code> is strictly greater than <code>low</code>.
 * This function is vectorized and will check each element of
 * <code>y</code> against each element of <code>low</code>.
 * @tparam T_y Type of y
 * @tparam T_low Type of lower bound
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw <code>domain_error</code> if y is not greater than low or
 *   if any element of y or low is NaN.
 */
template <typename T_y, typename T_low, require_container_t<T_y>* = nullptr>
inline void check_greater(const char* function, const char* name, const T_y& y,
                          const T_low& low) {
  const auto& low_ref = to_ref(value_of(low));
  scalar_seq_view<decltype(low_ref)> low_vec(low_ref);
  const auto& y_ref = to_ref(value_of(y));
  for (size_t n = 0; n < stan::math::size(y_ref); n++) {
    if (!(stan::get(y_ref, n) > low_vec[n])) {
      [&]() STAN_COLD_PATH {
        std::stringstream msg;
        msg << ", but must be greater than ";
        msg << low_vec[n];
        std::string msg_str(msg.str());
        throw_domain_error_vec(function, name, y_ref, n, "is ",
                               msg_str.c_str());
      }();
    }
  }
}

/**
 * Check if <code>y</code> is strictly greater than <code>low</code>.
 * This function is vectorized and will check each element of
 * <code>y</code> against each element of <code>low</code>.
 * @tparam T_y Type of y
 * @tparam T_low Type of lower bound
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw <code>domain_error</code> if y is not greater than low or
 *   if any element of y or low is NaN.
 */
template <typename T_y, typename T_low, require_not_container_t<T_y>* = nullptr>
inline void check_greater(const char* function, const char* name, const T_y& y,
                          const T_low& low) {
  const auto& low_ref = to_ref(value_of(low));
  scalar_seq_view<decltype(low_ref)> low_vec(low_ref);
  for (size_t n = 0; n < stan::math::size(low); n++) {
    if (!(y > low_vec[n])) {
      [&]() STAN_COLD_PATH {
        std::stringstream msg;
        msg << ", but must be greater than ";
        msg << low_vec[n];
        std::string msg_str(msg.str());
        throw_domain_error(function, name, y, "is ", msg_str.c_str());
      }();
    }
  }
}

template <typename T_y, typename T_low,
          require_any_var_matrix_t<T_y, T_low>* = nullptr>
inline void check_greater(const char* function, const char* name, const T_y& y,
                          const T_low& low) {
  check_greater(function, name, value_of(y), value_of(low));
}

}  // namespace math
}  // namespace stan
#endif
