#ifndef STAN_MATH_PRIM_ERR_CHECK_POSITIVE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_POSITIVE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/elementwise_check.hpp>
#include <stan/math/prim/err/invalid_argument.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <type_traits>
#include <string>

namespace stan {
namespace math {

/**
 * Check if <code>y</code> is positive.
 * This function is vectorized and will check each element of
 * <code>y</code>.
 * @tparam T_y Type of y
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param y Variable to check
 * @throw <code>domain_error</code> if y is negative or zero or
 *   if any element of y is NaN
 */
template <typename T_y>
inline void check_positive(const char* function, const char* name,
                           const T_y& y) {
  elementwise_check([](double x) { return x > 0; }, function, name, y,
                    "positive");
}

/**
 * Check if <code>size</code> is positive.
 * @param function Function name (for error messages)
 * @param name Variable name (for error messages)
 * @param expr Expression for the dimension size (for error messages)
 * @param size Size value to check
 * @throw <code>std::invalid_argument</code> if <code>size</code> is
 *   zero or negative.
 */
inline void check_positive(const char* function, const char* name,
                           const char* expr, int size) {
  if (size <= 0) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << "; dimension size expression = " << expr;
      std::string msg_str(msg.str());
      invalid_argument(function, name, size,
                       "must have a positive size, but is ", msg_str.c_str());
    }();
  }
}

}  // namespace math
}  // namespace stan
#endif
