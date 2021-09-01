#ifndef STAN_MATH_REV_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_REV_FUN_VALUE_OF_REC_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/value_of_rec.hpp>

namespace stan {
namespace math {

/**
 * Return the value of the specified variable.
 *
 * @param v Variable.
 * @return Value of variable.
 */
template <typename T, require_var_t<T>* = nullptr>
inline auto value_of_rec(const T& v) {
  return v.vi_->val_;
}

}  // namespace math
}  // namespace stan
#endif
