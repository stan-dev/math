#ifndef STAN_MATH_VERSION_HPP
#define STAN_MATH_VERSION_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <boost/version.hpp>
#include <sundials/sundials_config.h>
#include <string>

#if __has_include(<tbb/tbb_stddef.h>)
#include <tbb/tbb_stddef.h>
#else
#include <tbb/version.h>
#endif

#if EIGEN_VERSION_AT_LEAST(3, 4, 0)
#error "Stan Math is not yet compatible with Eigen 3.4"
#endif

// Hypergeometric functions were not introduced until Boost 1.72
#if BOOST_VERSION < 107200
#error "Boost 1.72 is the minimum compatible version for use with Stan Math"
#endif

#if SUNDIALS_VERSION_MAJOR < 6
#define SUNDIALS_INTERFACE_OLD
#endif

#if TBB_VERSION_MAJOR >= 2020
#define TBB_INTERFACE_NEW
#endif

#ifndef STAN_STRING_EXPAND
#define STAN_STRING_EXPAND(s) #s
#endif

#ifndef STAN_STRING
#define STAN_STRING(s) STAN_STRING_EXPAND(s)
#endif

#define STAN_MATH_MAJOR 4
#define STAN_MATH_MINOR 4
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
