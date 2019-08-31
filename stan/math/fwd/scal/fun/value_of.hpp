#ifndef STAN_MATH_FWD_SCAL_FUN_VALUE_OF_HPP
#define STAN_MATH_FWD_SCAL_FUN_VALUE_OF_HPP

#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {

/**
 * Return partial type for an fvar.
 * @tparam T The type of fvar.
 * @param x The fvar.
 * @return inner partial type of fvar.
 */
template <typename T, require_fvar<T>...>
inline auto&& value_of(T&& x) {
  return std::forward<T>(x).val_;
}

}  // namespace math
}  // namespace stan
#endif
