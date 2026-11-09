#ifndef STAN_MATH_PRIM_FUN_AS_VALUE_COLUMN_ARRAY_OR_SCALAR_REF_HPP
#define STAN_MATH_PRIM_FUN_AS_VALUE_COLUMN_ARRAY_OR_SCALAR_REF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <utility>

namespace stan {
namespace math {

/**
 * Return the values of a scalar or a vector as a scalar or an Eigen column
 * array. Use it in place of `to_ref(as_value_column_array_or_scalar(a))` to
 * avoid a copy of the values of an autodiff argument.
 *
 * For an argument without autodiff types the result is that of
 * `to_ref(as_value_column_array_or_scalar(a))`: an expression is evaluated
 * once, and a plain object is referenced. For an autodiff argument the result
 * is a view that reads the values from the autodiff variables at each access,
 * so the values are not copied.
 *
 * @tparam T a scalar, an Eigen column or row vector (or an expression of
 * one), a `std::vector` of scalars, or a `var_value` of an Eigen vector
 * @param a argument; an lvalue container must outlive the result
 * @return the values of `a`
 */
template <typename T>
inline decltype(auto) as_value_column_array_or_scalar_ref(T&& a) {
  return to_ref_if<!is_autodiff_v<T>>(
      as_value_column_array_or_scalar(std::forward<T>(a)));
}

}  // namespace math
}  // namespace stan

#endif
