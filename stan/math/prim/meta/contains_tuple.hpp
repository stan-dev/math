#ifndef STAN_MATH_PRIM_META_CONTAINS_TUPLE_HPP
#define STAN_MATH_PRIM_META_CONTAINS_TUPLE_HPP
#include <stan/math/prim/meta/is_vector.hpp>
#include <stan/math/prim/meta/is_tuple.hpp>

namespace stan {
namespace internal {
template <typename T>
struct contains_tuple_impl {
  static constexpr bool value = stan::math::is_tuple_v<T>;
};

template <typename T, typename... VecArgs>
struct contains_tuple_impl<std::vector<T, VecArgs...>> {
  static constexpr bool value = contains_tuple_impl<T>::value;
};
}  // namespace internal

/**
 * Check if the type is a tuple or contains a tuple
 * @tparam T type to check
 */
template <typename T>
struct contains_tuple {
  static constexpr bool value
      = internal::contains_tuple_impl<std::decay_t<T>>::value;
};
template <typename T>
inline constexpr bool contains_tuple_v = contains_tuple<std::decay_t<T>>::value;

/**
 * Check if the type is a std::vector containing a tuple
 */
template <typename T>
inline constexpr bool is_std_vector_containing_tuple_v
    = is_std_vector_v<std::decay_t<T>>&& contains_tuple_v<std::decay_t<T>>;

}  // namespace stan

#endif
