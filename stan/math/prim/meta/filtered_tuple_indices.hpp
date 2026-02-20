#ifndef STAN_MATH_PRIM_META_FILTERED_TUPLE_INDICES_HPP
#define STAN_MATH_PRIM_META_FILTERED_TUPLE_INDICES_HPP

#include <stan/math/prim/meta/filter_index_sequence.hpp>
#include <tuple>
#include <type_traits>
#include <utility>

namespace stan {
namespace math {

/**
 * Type trait returning the number of elements in an integer sequence.
 *
 * @tparam Sequence `std::integer_sequence` type
 * @return Inherits `std::integral_constant<std::size_t, N>` where `N` is the
 * sequence length.
 * @throw None.
 */
template <typename Sequence>
struct index_sequence_size;

/**
 * Specialization extracting the size of a concrete integer sequence.
 *
 * @tparam T integral value type used by the sequence
 * @tparam Is sequence values
 * @return Inherits `std::integral_constant<std::size_t, sizeof...(Is)>`.
 * @throw None.
 */
template <typename T, T... Is>
struct index_sequence_size<std::integer_sequence<T, Is...>>
    : std::integral_constant<std::size_t, sizeof...(Is)> {};

/**
 * Produces a compact index sequence of top-level tuple elements whose
 * types satisfy `Predicate`.
 *
 * @tparam Predicate unary predicate type trait
 * @tparam Tuple tuple type to filter
 * @return Provides nested `type` alias with filtered indices.
 * @throw None.
 */
template <template <class...> class Predicate, typename Tuple>
struct filtered_tuple_indices;

/**
 * Tuple specialization that filters top-level tuple indices using `Predicate`.
 *
 * @tparam Predicate unary predicate type trait
 * @tparam Args tuple element types
 * @return Provides nested `type` alias equal to the filtered index sequence.
 * @throw None.
 */
template <template <class...> class Predicate, typename... Args>
struct filtered_tuple_indices<Predicate, std::tuple<Args...>> {
  using type = filter_index_sequence_t<Predicate, std::index_sequence_for<Args...>,
                                       std::decay_t<Args>...>;
};

/**
 * Alias for the compact sequence of tuple indices whose element types satisfy
 * `Predicate`.
 *
 * @tparam Predicate unary predicate type trait
 * @tparam Tuple tuple type to filter
 * @return Filtered tuple index-sequence type.
 * @throw None.
 */
template <template <class...> class Predicate, typename Tuple>
using filtered_tuple_indices_t
    = typename filtered_tuple_indices<Predicate, std::decay_t<Tuple>>::type;

/**
 * Type trait returning the number of filtered tuple indices.
 *
 * @tparam Predicate unary predicate type trait
 * @tparam Tuple tuple type to filter
 * @return `std::integral_constant<std::size_t, N>` where `N` is the filtered
 * arity.
 * @throw None.
 */
template <template <class...> class Predicate, typename Tuple>
using filtered_tuple_size
    = index_sequence_size<filtered_tuple_indices_t<Predicate, Tuple>>;

/**
 * Convenience variable template for `filtered_tuple_size`.
 *
 * @tparam Predicate unary predicate type trait
 * @tparam Tuple tuple type to filter
 * @return Number of top-level tuple elements satisfying `Predicate`.
 * @throw None.
 */
template <template <class...> class Predicate, typename Tuple>
inline constexpr std::size_t filtered_tuple_size_v
    = filtered_tuple_size<Predicate, Tuple>::value;

}  // namespace math
}  // namespace stan

#endif
