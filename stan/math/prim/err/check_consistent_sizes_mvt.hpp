#ifndef STAN_MATH_PRIM_ERR_CHECK_CONSISTENT_SIZES_MVT_HPP
#define STAN_MATH_PRIM_ERR_CHECK_CONSISTENT_SIZES_MVT_HPP

#include <stan/math/prim/err/invalid_argument.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <algorithm>
#include <sstream>
#include <string>

namespace stan {
namespace math {

/** Trivial no input case, this function is a no-op. */
inline void check_consistent_sizes_mvt(const char*) { return; }

/**
 * Base case of recursion, this function is a no-op.
 * @tparam T1 type of first input
 **/
template <typename T1>
inline void check_consistent_sizes_mvt(const char*, const char*, const T1&) {
  return;
}

/**
 * Check that the provided multivariate inputs are of consistent size with each
 * other. Two multivariate inputs are of consistent size if both are
 * `std::vector`s of the same size, or if at least one is a not an
 * `std::vector`.
 *
 * E.g.: check_consistent_sizes_mvt("some_function", "x1", x1, "x2", x2, etc.).
 *
 * @tparam T1 type of first input
 * @tparam T2 type of second input
 * @tparam Ts type of other inputs
 * @param function function name (for error messages)
 * @param name1 name of variable corresponding to first input
 * @param x1 first input
 * @param name2 name of variable corresponding to second input
 * @param x2 second input
 * @param names_and_xs more inputs
 * @throw `invalid_argument` if sizes are inconsistent
 */
template <typename T1, typename T2, typename... Ts>
inline void check_consistent_sizes_mvt(const char* function, const char* name1,
                                       const T1& x1, const char* name2,
                                       const T2& x2,
                                       const Ts&... names_and_xs) {
  if (!is_std_vector<T1>::value && is_std_vector<T2>::value) {
    check_consistent_sizes_mvt(function, name2, x2, name1, x1, names_and_xs...);
  } else if (!is_std_vector<T2>::value) {
    check_consistent_sizes_mvt(function, name1, x1, names_and_xs...);
  } else if (stan::math::size(x1) == stan::math::size(x2)) {
    check_consistent_sizes_mvt(function, name1, x1, names_and_xs...);
  } else {
    [&]() STAN_COLD_PATH {
      size_t size_x1 = stan::math::size(x1);
      size_t size_x2 = stan::math::size(x2);
      std::stringstream msg;
      msg << ", but " << name2 << " has size " << size_x2
          << "; and they must be the same size.";
      std::string msg_str(msg.str());
      invalid_argument(function, name1, size_x1,
                       "has size = ", msg_str.c_str());
    }();
  }
}

}  // namespace math
}  // namespace stan
#endif
