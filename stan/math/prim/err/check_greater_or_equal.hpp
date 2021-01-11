#ifndef STAN_MATH_PRIM_ERR_CHECK_GREATER_OR_EQUAL_HPP
#define STAN_MATH_PRIM_ERR_CHECK_GREATER_OR_EQUAL_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/throw_domain_error_vec.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/scalar_seq_view.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <string>

namespace stan {
namespace math {

namespace internal {
template <typename T_y, typename T_low, bool is_vec>
struct greater_or_equal {
  static void check(const char* function, const char* name, const T_y& y,
                    const T_low& low) {
    scalar_seq_view<T_low> low_vec(low);
    for (size_t n = 0; n < stan::math::size(low); n++) {
      if (!(y >= low_vec[n])) {
        [&]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be greater than or equal to ";
          msg << low_vec[n];
          std::string msg_str(msg.str());
          throw_domain_error(function, name, y, "is ", msg_str.c_str());
        }();
      }
    }
  }
};

template <typename T_y, typename T_low>
struct greater_or_equal<T_y, T_low, true> {
  static void check(const char* function, const char* name, const T_y& y,
                    const T_low& low) {
    scalar_seq_view<T_low> low_vec(low);
    const auto& y_ref = to_ref(y);
    for (size_t n = 0; n < stan::math::size(y_ref); n++) {
      if (!(stan::get(y_ref, n) >= low_vec[n])) {
        [&]() STAN_COLD_PATH {
          std::stringstream msg;
          msg << ", but must be greater than or equal to ";
          msg << low_vec[n];
          std::string msg_str(msg.str());
          throw_domain_error_vec(function, name, y_ref, n, "is ",
                                 msg_str.c_str());
        }();
      }
    }
  }
};
}  // namespace internal

/**
 * Check if <code>y</code> is greater or equal than <code>low</code>.
 * This function is vectorized and will check each element of
 * <code>y</code> against each element of <code>low</code>.
 * @tparam T_y Type of y
 * @tparam T_low Type of lower bound
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @param low Lower bound
 * @throw <code>domain_error</code> if y is not greater or equal to low or
 *   if any element of y or low is NaN.
 */
template <typename T_y, typename T_low>
inline void check_greater_or_equal(const char* function, const char* name,
                                   const T_y& y, const T_low& low) {
  internal::greater_or_equal<T_y, T_low, is_vector_like<T_y>::value>::check(
      function, name, y, low);
}
}  // namespace math
}  // namespace stan
#endif
