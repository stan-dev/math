#ifndef STAN_MATH_PRIM_META_MODIFY_NESTED_VALUE_TYPE_HPP
#define STAN_MATH_PRIM_META_MODIFY_NESTED_VALUE_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_container.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Template metaprogram for modifying the nested types of containers.
 * This is the base case.
 *
 * @tparam TypeModifier Templated struct or alias that takes and returns
 *                        a single type
 * @tparam S Input type to modified
 */
template <template <class...> class TypeModifier, typename S,
          typename Enable = void>
struct modify_nested_value_type {
  using type = TypeModifier<S>;
};

/**
 * Template metaprogram for modifying the nested types of containers.
 *
 * Overload for std::vectors of container types
 *
 * @tparam TypeModifier Templated struct or alias that takes and returns
 *                        a single type
 * @tparam S Input type to modified
 */
template <template <class...> class TypeModifier, typename S>
struct modify_nested_value_type<TypeModifier, std::vector<S>,
                                require_container_t<S>> {
  using type = std::vector<TypeModifier<S>>;
};

/**
 * Template metaprogram for modifying the nested types of containers.
 *
 * Overload for tuples.
 *
 * @tparam TypeModifier Templated struct or alias that takes and returns
 *                        a single type
 * @tparam S Input type to modified
 */
template <template <class...> class TypeModifier, typename... S>
struct modify_nested_value_type<TypeModifier, std::tuple<S...>> {
  using type = std::tuple<TypeModifier<S>...>;
};

template <template <class...> class TypeModifier, typename S>
using modify_nested_value_type_t =
    typename modify_nested_value_type<TypeModifier, std::decay_t<S>>::type;

}  // namespace math
}  // namespace stan
#endif
