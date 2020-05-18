#ifndef STAN_MATH_REV_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_REV_FUN_VALUE_OF_REC_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun.hpp>
#include <stan/math/rev/fun/value_of.hpp>

namespace stan {
namespace math {

/**
 * Return the value of the specified variable.
 *
 * @param v Variable.
 * @return Value of variable.
 */
inline auto value_of_rec(const var& v) { return v.vi_->val_; }

template <typename Vec, require_std_vector_vt<is_var, Vec>* = nullptr>
inline auto value_of_rec(Vec&& x) {
  return value_of(std::forward<Vec>(x));
}

}  // namespace math
}  // namespace stan
#endif
