#ifndef STAN_MATH_PRIM_META_TUPLE_ELEMENT_HPP
#define STAN_MATH_PRIM_META_TUPLE_ELEMENT_HPP

#include <stan/math/prim/meta/is_tuple.hpp>
#include <stan/math/prim/meta/tuple_size.hpp>

#include <type_traits>
#include <cstddef>
#include <tuple>

namespace stan {

/**
 * Equivalent to std::tuple_element but returns void if N is out of range
 * to avoid a static assertion failure in libstdc++.
 * @tparam N index of the element to retrieve
 * @tparam T type to retrieve the element from
 * @ingroup type_trait
 */
template <std::size_t N, typename T, typename = void>
struct tuple_element {
  using type = void;
};

template <std::size_t N, typename T>
struct tuple_element<
    N, T,
    std::enable_if_t<stan::is_tuple_v<T> && (N < stan::tuple_size_v<T>)>> {
  using type = std::tuple_element_t<N, std::decay_t<T>>;
};

template <std::size_t N, typename T>
using tuple_element_t = typename stan::tuple_element<N, T>::type;
}  // namespace stan

#endif
