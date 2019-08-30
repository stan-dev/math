#ifndef STAN_MATH_REV_SCAL_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_REV_SCAL_FUN_VALUE_OF_REC_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/scal/fun/value_of_rec.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Return the value of the specified variable.
 *
 * @param x Variable.
 * @return Value of variable.
 */
template <typename T, enable_if_var<std::decay_t<T>>* = nullptr>
inline auto&& value_of_rec(T&& x) {
  return x.vi_->val_;
}

}  // namespace math
}  // namespace stan
#endif
