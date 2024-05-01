#ifndef STAN_MATH_PRIM_META_MODIFY_NESTED_VALUE_TYPE_HPP
#define STAN_MATH_PRIM_META_MODIFY_NESTED_VALUE_TYPE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_container.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Template metaprogram to calculate a type for converting a
 * convertible type.  This is the base case.
 *
 * @tparam T result scalar type.
 * @tparam S input type
 */
template <template <class...> class TypeModifier, typename S,
          typename Enable = void>
struct modify_nested_value_type {
  /**
   * The promoted type.
   */
  using type = TypeModifier<S>;
};

/**
 * Template metaprogram to calculate a type for a container whose
 * underlying scalar is converted from the second template
 * parameter type to the first.
 *
 * @tparam T result scalar type.
 * @tparam S input type
 */
template <template <class...> class TypeModifier, typename S>
struct modify_nested_value_type<TypeModifier, std::vector<S>,
                                require_container_t<S>> {
  /**
   * The promoted type.
   */
  using type =
    std::vector<typename modify_nested_value_type<TypeModifier, S>::type>;
};

template <template <class...> class TypeModifier, typename... S>
struct modify_nested_value_type<TypeModifier, std::tuple<S...>> {
  /**
   * The promoted type.
   */
  using type =
    std::tuple<typename modify_nested_value_type<TypeModifier, S>::type...>;
};

template <template <class...> class TypeModifier, typename S>
using modify_nested_value_type_t =
    typename modify_nested_value_type<TypeModifier, std::decay_t<S>>::type;

}  // namespace math
}  // namespace stan
#endif
