#ifndef STAN_MATH_PRIM_ERR_CHECK_FINITE_HPP
#define STAN_MATH_PRIM_ERR_CHECK_FINITE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/elementwise_check.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/get.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if all values in `y` are finite. `y` can be a
 *scalar, `std::vector` or Eigen type.
 *
 * @tparam T_y type of `y`
 *
 * @param function name of function (for error messages)
 * @param name variable name (for error messages)
 * @param y scalar or container to test
 * @return <code>true</code> if all values are finite
 **/
template <typename T_y>
inline void check_finite(const char* function, const char* name, const T_y& y) {
  elementwise_check([](double x) { return std::isfinite(x); }, function, name,
                    y, "finite");
}

}  // namespace math
}  // namespace stan

#endif
