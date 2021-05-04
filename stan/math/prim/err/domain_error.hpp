#ifndef STAN_MATH_PRIM_ERR_DOMAIN_ERROR_HPP
#define STAN_MATH_PRIM_ERR_DOMAIN_ERROR_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/throw_domain_error.hpp>

namespace stan {
namespace math {

/**
 * @deprecated use <code>throw_domain_error</code>
 */
template <typename T>
inline void domain_error(const char* function, const char* name, const T& y,
                         const char* msg1, const char* msg2) {
  throw_domain_error(function, name, y, msg1, msg2);
}

/**
 * @deprecated use <code>throw_domain_error</code>
 */
template <typename T>
inline void domain_error(const char* function, const char* name, const T& y,
                         const char* msg1) {
  throw_domain_error(function, name, y, msg1);
}

}  // namespace math
}  // namespace stan
#endif
