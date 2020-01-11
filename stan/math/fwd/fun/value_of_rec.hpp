#ifndef STAN_MATH_FWD_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_FWD_FUN_VALUE_OF_REC_HPP

#include <stan/math/prim/fun/value_of_rec.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>

namespace stan {
namespace math {

/**
 * Return the value of the specified variable.
 *
 * @tparam T inner type of the fvar, must implement value_of_rec
 * @param v Variable.
 * @return Value of variable.
 */

template <typename T>
inline double value_of_rec(const fvar<T>& v) {
  return value_of_rec(v.val_);
}

}  // namespace math
}  // namespace stan
#endif
