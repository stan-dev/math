#ifndef STAN_MATH_PRIM_ERR_CHECK_CONSISTENT_SIZE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_CONSISTENT_SIZE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/invalid_argument.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <sstream>
#include <string>

namespace stan {
namespace math {

/**
 * Check if `x` is consistent with size `expected_size`. `x` is consistent with
 * size `expected_size` if `x` is a vector of size `expected_size`, or a non
 * vector.
 * @tparam T type of value
 * @param function function name (for error messages)
 * @param name variable name (for error messages)
 * @param x variable to check for consistent size
 * @param expected_size expected size if `x` is a vector
 * @throw `invalid_argument` if the size is inconsistent
 */
template <typename T>
inline void check_consistent_size(const char* function, const char* name,
                                  const T& x, size_t expected_size) {
  if (!(!is_vector<T>::value
        || (is_vector<T>::value && expected_size == stan::math::size(x)))) {
    [&]() STAN_COLD_PATH {
      std::stringstream msg;
      msg << ", expecting dimension = " << expected_size
          << "; a function was called with arguments of different "
          << "scalar, array, vector, or matrix types, and they were not "
          << "consistently sized;  all arguments must be scalars or "
          << "multidimensional values of the same shape.";
      std::string msg_str(msg.str());

      invalid_argument(function, name, stan::math::size(x),
                       "has dimension = ", msg_str.c_str());
    }();
  }
}

}  // namespace math
}  // namespace stan
#endif
