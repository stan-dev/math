#ifndef STAN_MATH_PRIM_MAT_FUN_VALUE_OF_HPP
#define STAN_MATH_PRIM_MAT_FUN_VALUE_OF_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/scal/fun/value_of.hpp>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {

/**
 * Return the specified argument.
 *
 * <p>See <code>value_of(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @param x Specified matrix.
 * @return Specified matrix.
 */
template <typename T, enable_if_eigen<T>...,
          enable_if_floating_point<scalar_type_decay_t<T>>...>
inline auto&& value_of(T&& x) {
  return std::forward<T>(x);
}

template <typename T, enable_if_eigen<T>...,
          enable_if_arithmetic<scalar_type_decay_t<T>>..., enable_if_not_floating_point<scalar_type_decay_t<T>>...>
inline auto value_of(T&& x) {
  return (x.template cast<double>()).eval();
}

}  // namespace math
}  // namespace stan

#endif
