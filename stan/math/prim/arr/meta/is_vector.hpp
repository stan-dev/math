#ifndef STAN_MATH_PRIM_ARR_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_ARR_META_IS_VECTOR_HPP

#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <vector>
#include <type_traits>
namespace stan {

namespace internal {
  // FIXME(Steve): this is a stupid name
  template <typename T>
  using inherit_vector = std::vector<typename std::decay_t<T>::value_type, typename std::decay_t<T>::allocator_type>;
}

/**
 * Metaprogram structure to determine if type is a standard vector
 *
 * <p>This base class should be specialized for structured types.</p>
 *
 * @tparam T Type of object.
 */
template <typename T, typename = void>
struct is_std_vector : std::false_type {};

/**
 * Specialization of is_std_vector for standard vectors.
 *
 * @tparam T Type of object.
 */
template <typename T>
struct is_std_vector<T,
        typename std::enable_if_t<std::is_same<typename std::decay_t<T>, internal::inherit_vector<T>>::value>> : std::true_type {
  typedef std::decay_t<T> type;
};

/**
 * Specialization of is_vector for standard vectors.
 *
 * @tparam T Type of object.
 */
template<typename T>
struct is_vector<T,
        typename std::enable_if_t<std::is_same<typename std::decay_t<T>, internal::inherit_vector<T>>::value>> : is_std_vector<T> {
  typedef std::decay_t<T> type;
};



}  // namespace stan
#endif
