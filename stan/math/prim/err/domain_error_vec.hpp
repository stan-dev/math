#ifndef STAN_MATH_PRIM_ERR_DOMAIN_ERROR_VEC_HPP
#define STAN_MATH_PRIM_ERR_DOMAIN_ERROR_VEC_HPP

#include <stan/math/prim/err/throw_domain_error_vec.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>throw_domain_error_vec</code>
 */
template <typename T>
inline void domain_error_vec(const char* function, const char* name, const T& y,
                             size_t i, const char* msg1, const char* msg2) {
  throw_domain_error_vec(function, name, y, i, msg1, msg2);
}

/**
 * @deprecated use <code>throw_domain_error_vec</code>
 */
template <typename T>
inline void domain_error_vec(const char* function, const char* name, const T& y,
                             size_t i, const char* msg1) {
  throw_domain_error_vec(function, name, y, i, msg1, "");
}

}  // namespace math
}  // namespace stan
#endif
