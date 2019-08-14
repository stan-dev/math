#ifndef STAN_MATH_PRIM_SCAL_META_ENABLE_IF_SCALAR_ARITHMETIC_HPP
#define STAN_MATH_PRIM_SCAL_META_ENABLE_IF_SCALAR_ARITHMETIC_HPP

#include <stan/math/prim/scal/meta/conjunction.hpp>
#include <stan/math/prim/scal/meta/disjunction.hpp>
#include <stan/math/prim/scal/meta/scalar_type.hpp>
#include <type_traits>

namespace stan {

template <typename T>
using enable_if_scalar_arithmetic = std::enable_if_t<
    std::is_arithmetic<scalar_type_t<std::decay_t<T>>>::value>;

template <typename T>
using enable_if_not_scalar_arithmetic = std::enable_if_t<
    !std::is_arithmetic<scalar_type_t<std::decay_t<T>>>::value>;

template <typename... Types>
using enable_if_all_scalar_arithmetic = std::enable_if_t<math::conjunction<
    std::is_arithmetic<scalar_type_t<std::decay_t<Types>>>...>::value>;

template <typename... Types>
using enable_if_any_scalar_arithmetic = std::enable_if_t<math::disjunction<
    std::is_arithmetic<scalar_type_t<std::decay_t<Types>>>...>::value>;

template <typename... Types>
using enable_if_all_not_scalar_arithmetic = std::enable_if_t<!math::conjunction<
    std::is_arithmetic<scalar_type_t<std::decay_t<Types>>>...>::value>;

template <typename... Types>
using enable_if_any_not_scalar_arithmetic = std::enable_if_t<!math::disjunction<
    std::is_arithmetic<scalar_type_t<std::decay_t<Types>>>...>::value>;

}  // namespace stan
#endif
