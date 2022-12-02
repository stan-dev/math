#ifndef STAN_MATH_VERSION_HPP
#define STAN_MATH_VERSION_HPP

#include <string>

#ifndef STAN_STRING_EXPAND
#define STAN_STRING_EXPAND(s) #s
#endif

#ifndef STAN_STRING
#define STAN_STRING(s) STAN_STRING_EXPAND(s)
#endif

#define STAN_MATH_MAJOR 4
#define STAN_MATH_MINOR 5
#define STAN_MATH_PATCH 0

namespace stan {
namespace math {

/** Major version number for Stan math library. */
const std::string MAJOR_VERSION = STAN_STRING(STAN_MATH_MAJOR);

/** Minor version number for Stan math library. */
const std::string MINOR_VERSION = STAN_STRING(STAN_MATH_MINOR);

/** Patch version for Stan math library. */
const std::string PATCH_VERSION = STAN_STRING(STAN_MATH_PATCH);

}  // namespace math
}  // namespace stan

#endif
