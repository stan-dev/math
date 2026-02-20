#ifndef STAN_MATH_PRIM_META_FILTERED_TUPLE_INDICES_HPP
#define STAN_MATH_PRIM_META_FILTERED_TUPLE_INDICES_HPP

#include <stan/math/prim/meta/filter_index_sequence.hpp>
#include <tuple>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {

template <typename Sequence>
struct index_sequence_size;

template <typename T, T... Is>
struct index_sequence_size<std::integer_sequence<T, Is...>>
    : std::integral_constant<std::size_t, sizeof...(Is)> {};

/**
 * Produces a compact index sequence of top-level tuple elements whose
 * types satisfy `Predicate`.
 */
template <template <class...> class Predicate, typename Tuple>
struct filtered_tuple_indices;

template <template <class...> class Predicate, typename... Args>
struct filtered_tuple_indices<Predicate, std::tuple<Args...>> {
  using type = filter_index_sequence_t<Predicate, std::index_sequence_for<Args...>,
                                       std::decay_t<Args>...>;
};

template <template <class...> class Predicate, typename Tuple>
using filtered_tuple_indices_t
    = typename filtered_tuple_indices<Predicate, std::decay_t<Tuple>>::type;

template <template <class...> class Predicate, typename Tuple>
using filtered_tuple_size
    = index_sequence_size<filtered_tuple_indices_t<Predicate, Tuple>>;

template <template <class...> class Predicate, typename Tuple>
inline constexpr std::size_t filtered_tuple_size_v
    = filtered_tuple_size<Predicate, Tuple>::value;

}  // namespace math
}  // namespace stan

#endif
