#ifndef STAN_MATH_PRIM_SCAL_META_IS_VECTOR_LIKE_HPP
#define STAN_MATH_PRIM_SCAL_META_IS_VECTOR_LIKE_HPP

#include <stan/math/prim/scal/meta/is_vector.hpp>
#include <type_traits>

namespace stan {

/**
 * Template metaprogram indicates whether a type is vector_like.
 *
 * A type is vector_like if an instance can be accessed like a
 * vector, i.e. square brackets.
 *
 * Access is_vector_like::value for the result.
 *
 * Default behavior is to use the is_vector template metaprogram.
 *
 * @tparam T Type to test
 */
template <typename T>
struct is_vector_like
    : std::integral_constant<bool, stan::is_vector<T>::value
                                       || is_eigen<T>::value> {};

/**
 * Template metaprogram indicates whether a type is vector_like.
 *
 * A type is vector_like if an instance can be accessed like a
 * vector, i.e. square brackets.
 *
 * A C++ array of T is vector_like.
 *
 * @tparam T Type to test
 */
template <typename T>
struct is_vector_like<T*> : std::true_type {};

/**
 * Template metaprogram indicates whether a type is vector_like.
 *
 * A type is vector_like if an instance can be accessed like a
 * vector, i.e. square brackets.
 *
 * Access is_vector_like::value for the result.
 *
 * This metaprogram removes the const qualifier.
 *
 * @tparam T Type to test
 */
template <typename T>
struct is_vector_like<const T>
    : std::integral_constant<bool, stan::is_vector_like<T>::value> {};

}  // namespace stan
#endif
