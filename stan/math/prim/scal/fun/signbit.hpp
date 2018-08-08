#ifndef STAN_MATH_PRIM_SCAL_FUN_SIGNBIT_HPP
#define STAN_MATH_PRIM_SCAL_FUN_SIGNBIT_HPP

#include <stan/math/prim/scal/fun/val.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Determines if the argument is negative
 *
 * Needed for libc++'s implementation of
 * trig on complex values
 *
 * @tparam T type of stan AD object
 * @param t stan AD object
 * @return true if t is negative
 */
template <class T>
inline auto signbit(T const& t) {
  using std::signbit;
  return signbit(val(t));
}

}  // namespace math
}  // namespace stan
#endif
