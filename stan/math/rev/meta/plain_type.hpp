#ifndef STAN_MATH_REV_META_PLAIN_TYPE_HPP
#define STAN_MATH_REV_META_PLAIN_TYPE_HPP

#include <stan/math/prim/meta/is_eigen.hpp>
#include <stan/math/rev/meta/is_var.hpp>
#include <type_traits>

namespace stan {

/**
 * Determines plain (non expression) type associated with \c T. For `var_value`
 * with an underlying eigen type this is a `var_value<plain_type_t<T>>`
 *
 * @tparam T type to determine plain type of
 */
template <typename T>
struct plain_type<
    T,
    require_t<stan::math::conjunction<is_var<T>, is_eigen<value_type_t<T>>>>> {
  using type = typename stan::math::var_value<plain_type_t<value_type_t<T>>>;
};

}  // namespace stan

#endif  // STAN_MATH_PRIM_META_PLAIN_TYPE_HPP
