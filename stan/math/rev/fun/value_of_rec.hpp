#ifndef STAN_MATH_REV_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_REV_FUN_VALUE_OF_REC_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Return the value of the specified variable.
 *
 * @param v Variable.
 * @return Value of variable.
 */
template <typename T>
inline auto& value_of_rec(const var_value<T>& v) {
  return v.vi_->val_;
}

}  // namespace math
}  // namespace stan
#endif
