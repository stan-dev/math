#ifndef STAN_MATH_REV_SCAL_META_IS_VAR_HPP
#define STAN_MATH_REV_SCAL_META_IS_VAR_HPP

#include <stan/math/prim/scal/meta/is_var.hpp>
#include <stan/math/rev/core.hpp>
#include <type_traits>

namespace stan {
/**
 * Defines a public enum named value and sets it to true(1)
 * when instantiated with the stan::math::var type.
 */
template <typename T>
struct is_var<T, std::enable_if_t<std::is_same<stan::math::var, std::decay_t<T>>::value>> : std::true_type {};

}  // namespace stan
#endif
