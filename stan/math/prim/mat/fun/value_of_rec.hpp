#ifndef STAN_MATH_PRIM_MAT_FUN_VALUE_OF_REC_HPP
#define STAN_MATH_PRIM_MAT_FUN_VALUE_OF_REC_HPP

#include <stan/math/prim/scal/fun/value_of_rec.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Return the specified argument.
 *
 * <p>See <code>value_of_rec(T)</code> for a polymorphic
 * implementation using static casts.
 *
 * <p>This inline pass-through no-op should be compiled away.
 *
 * @param x Specified matrix.
 * @return Specified matrix.
 */
 template <typename T, enable_if_eigen<T>* = nullptr, enable_if_arithmetic<scalar_type_t<T>>* = nullptr>
 inline const auto&& value_of_rec(T&& x) {
   return x;
 }
}  // namespace math
}  // namespace stan

#endif
