#ifndef STAN_MATH_PRIM_META_REQUIRE_HELPERS_HPP
#define STAN_MATH_PRIM_META_REQUIRE_HELPERS_HPP

#include <stan/math/prim/meta/conjunction.hpp>
#include <stan/math/prim/meta/disjunction.hpp>
#include <stan/math/prim/meta/bool_constant.hpp>
#include <stan/math/prim/meta/scalar_type.hpp>
#include <stan/math/prim/meta/value_type.hpp>

#include <type_traits>

namespace stan {

/**
 * If condition is true, template is enabled
 * @ingroup type_traits
 */
template <class Check>
using require_t = std::enable_if_t<Check::value>;

/**
 * If condition is false, template is disabled
 * @ingroup type_traits
 */
template <typename Check>
using require_not_t = std::enable_if_t<!Check::value>;

/**
 * If all conditions are true, template is enabled
 * Returns a type void if all conditions are true and otherwise fails.
 * @ingroup type_traits
 */
template <class... Checks>
using require_all_t = std::enable_if_t<math::conjunction<Checks...>::value>;

/**
 * If any condition is true, template is enabled.
 *
 * Returns a type void if any of the conditions are true and otherwise fails.
 * @ingroup type_traits
 */
template <class... Checks>
using require_any_t = std::enable_if_t<math::disjunction<Checks...>::value>;

/**
 * If all conditions are false, template is enabled.
 *
 * Returns a type void if all of the conditions are false.
 * @ingroup type_traits
 */
template <class... Checks>
using require_all_not_t
    = std::enable_if_t<!math::disjunction<Checks...>::value>;

/**
 * If any condition is false, template is enabled.
 *
 * Returns a type void if any of the conditions are false.
 * @ingroup type_traits
 */
template <class... Checks>
using require_any_not_t
    = std::enable_if_t<!math::conjunction<Checks...>::value>;

/**
 * Used as the base for checking whether a type is a container with
 * an underlying scalar type
 *
 * @tparam ContainerCheck Templated struct or alias that wraps a static constant
 * scalar called type. Used to check the container satisfies a particular type
 * check.
 * @tparam ValueCheck Templated struct or alias that returns a containers
 * inner type.
 * @tparam CheckType Templated struct or alias that wraps a static constant
 * scalar called type. Used to check the container's underlying type satisfies a
 * particular type check.
 *
 */
template <template <class...> class ContainerCheck,
          template <class...> class ValueCheck,
          template <class...> class TypeCheck, class... Check>
using container_type_check_base
    = bool_constant<math::conjunction<ContainerCheck<std::decay_t<Check>>...,
                                      TypeCheck<ValueCheck<Check>>...>::value>;

}  // namespace stan

#endif
