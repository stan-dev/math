#ifndef STAN_MATH_PRIM_META_IS_TUPLE_ELT_HPP
#define STAN_MATH_PRIM_META_IS_TUPLE_ELT_HPP

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
    N, T, std::enable_if_t<(N < std::tuple_size_v<std::decay_t<T>>)>> {
  using type = std::tuple_element_t<N, std::decay_t<T>>;
};

template <std::size_t N, typename T>
using tuple_element_t = typename tuple_element<N, T>::type;
}  // namespace stan

#endif
