#ifndef STAN_MATH_REV_CORE_FILTER_VAR_SCALAR_TYPES_HPP
#define STAN_MATH_REV_CORE_FILTER_VAR_SCALAR_TYPES_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/functor/filter_map.hpp>
#include <stan/math/rev/meta.hpp>
#include <type_traits>

namespace stan::math::internal {
/**
 * Filter a tuple and return a tuple with references to the types with a var
 * scalar type.
 * @tparam T Possibly a tuple, std::vector, Eigen type, or scalar
 * @param[in] t Input to filter
 * @return Filtered input with only var scalar types
 */
template <typename T>
inline constexpr decltype(auto) filter_var_scalar_types(T&& t) {
  return stan::math::filter_map<is_any_var_scalar>(
      [](auto&& arg) -> decltype(auto) {
        return std::forward<decltype(arg)>(arg);
      },
      std::forward<T>(t));
}

}  // namespace stan::math::internal
#endif
