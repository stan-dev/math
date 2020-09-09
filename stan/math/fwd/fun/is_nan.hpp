#ifndef STAN_MATH_FWD_FUN_IS_NAN_HPP
#define STAN_MATH_FWD_FUN_IS_NAN_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/is_nan.hpp>

namespace stan {
namespace math {

/**
 * Returns 1 if the input's value is NaN and 0 otherwise.
 *
 * Delegates to <code>is_nan</code>.
 *
 * @tparam T inner type of the fvar
 * @param x Value to test.
 * @return <code>1</code> if the value is NaN and <code>0</code> otherwise.
 */
template <typename T, require_fvar_t<T>* = nullptr>
inline bool is_nan(T&& x) {
  return is_nan(std::forward<decltype(x.val())>(x.val()));
}

}  // namespace math
}  // namespace stan
#endif
