#ifndef STAN_MATH_PRIM_SCAL_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_VECTOR_HPP

#include <type_traits>

namespace stan {

/**
 * Metaprogram structure to determine if type is a vector
 *
 * <p>This base class should be specialized for structured types.</p>
 *
 * @tparam T Type of object.
 */
template <typename T, typename = void>
struct is_vector : std::false_type {
  typedef std::decay_t<T> type;
};

/**
 * Metaprogram structure to determine if type is a standard vector
 *
 * <p>This base class should be specialized for structured types.</p>
 *
 * @tparam T Type of object.
 */
template <typename T, typename = void>
struct is_std_vector : std::false_type {};

}  // namespace stan
#endif
