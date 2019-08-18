#ifndef STAN_MATH_PRIM_SCAL_META_IS_VECTOR_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_VECTOR_HPP

#include <stan/math/prim/scal/meta/is_eigen.hpp>
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

// Check whether type is an eigen col vector
template <typename T, typename = void>
struct is_eigen_col_vector : std::false_type {};


// Check whether type is an eigen row vector
template <typename T, typename = void>
struct is_eigen_row_vector : std::false_type {};

// Checks whether decayed type is an eigen vector
template <typename T, typename = void>
struct is_eigen_vector : std::false_type {};


}  // namespace stan
#endif
