#ifndef STAN_MATH_FWD_FUN_VALUE_OF_HPP
#define STAN_MATH_FWD_FUN_VALUE_OF_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>

namespace stan {
namespace math {

/**
 * Return the value of the specified variable.
 *
 * @tparam T inner type of the fvar
 * @param v Variable.
 * @return Value of variable.
 */
template <typename T>
inline T value_of(const fvar<T>& v) {
  return v.val_;
}

}  // namespace math
}  // namespace stan
#endif
