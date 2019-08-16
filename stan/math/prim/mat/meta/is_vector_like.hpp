#ifndef STAN_MATH_PRIM_MAT_META_IS_VECTOR_LIKE_HPP
#define STAN_MATH_PRIM_MAT_META_IS_VECTOR_LIKE_HPP

#include <stan/math/prim/scal/meta/is_vector_like.hpp>
#include <stan/math/prim/mat/meta/is_eigen.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>

namespace stan {

/**
 * Template metaprogram indicates whether a type is vector_like.
 *
 * A type is vector_like if an instance can be accessed like a
 * vector, i.e. square brackets.
 *
 * Access is_vector_like::value for the result.
 *
 * @tparam T Type to test
 */
template <typename T>
struct is_vector_like<T, std::enable_if_t<is_eigen_decay<T>::value>>
    : std::true_type {};

}  // namespace stan
#endif
