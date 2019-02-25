#ifndef STAN_MATH_PRIM_META_IS_VECTOR_LIKE_HPP
#define STAN_MATH_PRIM_META_IS_VECTOR_LIKE_HPP

#include <stanh/prim/meta/is_vector_like.hpp>
#include <stanh/prim/fun/Eigen.hpp>

namespace stan {

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
struct is_vector_like<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > {
  enum { value = true };
};
}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_META_IS_VECTOR_LIKE_HPP
#define STAN_MATH_PRIM_META_IS_VECTOR_LIKE_HPP

#include <stanh/prim/meta/is_vector_like.hpp>
#include <stanh/prim/fun/Eigen.hpp>

namespace stan {

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
struct is_vector_like<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > {
  enum { value = true };
};
}  // namespace stan
#endif
