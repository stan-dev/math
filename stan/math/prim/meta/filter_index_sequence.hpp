#ifndef STAN_MATH_PRIM_META_FILTER_INDEX_SEQUENCE_HPP
#define STAN_MATH_PRIM_META_FILTER_INDEX_SEQUENCE_HPP

#include <utility>

namespace stan {
namespace math {

/**
 * Concatenate two integer sequences of the same value type.
 *
 * @tparam SeqT integral value type used by the sequences
 * @tparam Is first sequence indices
 * @tparam Js second sequence indices
 * @param[in] first first integer sequence
 * @param[in] second second integer sequence
 * @return concatenated sequence containing `Is...` followed by `Js...`
 * @throw None.
 */
template <typename SeqT, SeqT... Is, SeqT... Js>
inline constexpr std::integer_sequence<SeqT, Is..., Js...>
concat_integer_sequences([[maybe_unused]] std::integer_sequence<SeqT, Is...> first,
                         [[maybe_unused]] std::integer_sequence<SeqT, Js...> second) {
  return {};
}

/**
 * Identity overload for single-sequence concatenation.
 *
 * @tparam SeqT integral value type used by the sequence
 * @param[in] seq input integer sequence
 * @return `seq` unchanged
 * @throw None.
 */
template <typename SeqT>
inline constexpr std::integer_sequence<SeqT> concat_integer_sequences(
    [[maybe_unused]] std::integer_sequence<SeqT> seq) {
  return {};
}

/**
 * Concatenate three or more integer sequences.
 *
 * @tparam SeqT integral value type used by the sequences
 * @tparam Is indices in the first sequence
 * @tparam Js indices in the second sequence
 * @tparam Rest remaining sequence types
 * @param[in] first first integer sequence
 * @param[in] second second integer sequence
 * @param[in] rest remaining integer sequences
 * @return concatenated sequence preserving argument order
 * @throw None.
 */
template <typename SeqT, SeqT... Is, SeqT... Js, class... Rest>
inline constexpr auto concat_integer_sequences(
    [[maybe_unused]] std::integer_sequence<SeqT, Is...> first,
    [[maybe_unused]] std::integer_sequence<SeqT, Js...> second,
    [[maybe_unused]] Rest... rest) {
  return concat_integer_sequences(std::integer_sequence<SeqT, Is..., Js...>{},
                                  Rest{}...);
}

/**
 * Returns a singleton index sequence if `Predicate<Arg>::value` is true,
 * otherwise an empty sequence.
 *
 * @tparam Predicate unary predicate type trait
 * @tparam Arg argument type evaluated by `Predicate`
 * @tparam SeqT integral value type used by the sequence
 * @tparam I index candidate
 * @param[in] idx singleton input sequence containing `I`
 * @return `std::integer_sequence<SeqT, I>` if `Predicate<Arg>::value`,
 * otherwise `std::integer_sequence<SeqT>{}`
 * @throw None.
 */
template <template <class...> class Predicate, typename Arg, class SeqT,
          SeqT I>
inline constexpr auto filter_single_index_sequence(
    [[maybe_unused]] std::integer_sequence<SeqT, I> idx) {
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
 *
 * @tparam Predicate unary predicate type trait
 * @tparam Args argument type pack aligned with `Is...`
 * @tparam SeqT integral value type used by the sequence
 * @tparam Is index pack to filter
 * @param[in] idxs input index sequence
 * @return compact integer sequence of indices where
 * `Predicate<Args>::value` is true
 * @throw None. Fails compilation with `static_assert` if `sizeof...(Args) !=
 * sizeof...(Is)`.
 */
template <template <class...> class Predicate, typename... Args, class SeqT,
          SeqT... Is>
inline constexpr auto filter_index_sequence(
    [[maybe_unused]] std::integer_sequence<SeqT, Is...> idxs) {
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

/**
 * Alias returning the filtered sequence type produced by
 * `filter_index_sequence`.
 *
 * @tparam Predicate unary predicate type trait
 * @tparam Sequence input index sequence type
 * @tparam Args argument type pack aligned with `Sequence`
 */
template <template <class...> class Predicate, typename Sequence,
          typename... Args>
using filter_index_sequence_t
    = decltype(filter_index_sequence<Predicate, Args...>(Sequence{}));

}  // namespace math
}  // namespace stan

#endif
