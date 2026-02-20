#ifndef STAN_MATH_PRIM_META_FILTER_INDEX_SEQUENCE_HPP
#define STAN_MATH_PRIM_META_FILTER_INDEX_SEQUENCE_HPP

#include <utility>

namespace stan {
namespace math {

/**
 * Concatenates integer sequences.
 */
template <typename SeqT, SeqT... Is, SeqT... Js>
inline constexpr std::integer_sequence<SeqT, Is..., Js...>
concat_integer_sequences(std::integer_sequence<SeqT, Is...>,
                         std::integer_sequence<SeqT, Js...>) {
  return {};
}

template <typename SeqT>
inline constexpr std::integer_sequence<SeqT> concat_integer_sequences(
    std::integer_sequence<SeqT>) {
  return {};
}

template <typename SeqT, SeqT... Is, SeqT... Js, class... Rest>
inline constexpr auto concat_integer_sequences(
    std::integer_sequence<SeqT, Is...>, std::integer_sequence<SeqT, Js...>,
    Rest...) {
  return concat_integer_sequences(std::integer_sequence<SeqT, Is..., Js...>{},
                                  Rest{}...);
}

/**
 * Returns a singleton index sequence if `Predicate<Arg>::value` is true,
 * otherwise an empty sequence.
 */
template <template <class...> class Predicate, typename Arg, class SeqT,
          SeqT I>
inline constexpr auto filter_single_index_sequence(
    std::integer_sequence<SeqT, I>) {
  if constexpr (Predicate<Arg>::value) {
    return std::integer_sequence<SeqT, I>{};
  } else {
    return std::integer_sequence<SeqT>{};
  }
}

/**
 * Filters an index sequence by applying `Predicate` to matching `Args...`.
 *
 * `Args...` and the index sequence must have the same length.
 */
template <template <class...> class Predicate, typename... Args, class SeqT,
          SeqT... Is>
inline constexpr auto filter_index_sequence(std::integer_sequence<SeqT, Is...>) {
  static_assert(sizeof...(Args) == sizeof...(Is),
                "filter_index_sequence requires the same number of argument "
                "types and indices.");
  if constexpr (sizeof...(Is) > 0) {
    return concat_integer_sequences(
        filter_single_index_sequence<Predicate, Args>(
            std::integer_sequence<SeqT, Is>{})...);
  } else {
    return std::integer_sequence<SeqT>{};
  }
}

template <template <class...> class Predicate, typename Sequence,
          typename... Args>
using filter_index_sequence_t
    = decltype(filter_index_sequence<Predicate, Args...>(Sequence{}));

}  // namespace math
}  // namespace stan

#endif
