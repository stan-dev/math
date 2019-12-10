#ifndef STAN_MATH_PRIM_MAT_META_PLAIN_TYPE_HPP
#define STAN_MATH_PRIM_MAT_META_PLAIN_TYPE_HPP

#include <stan/math/prim/scal/meta/plain_type.hpp>
#include <stan/math/prim/scal/meta/require_generics.hpp>
#include <type_traits>

namespace stan {

/**
 * Determines plain (non expression) type associated with \c T. For \c Eigen
 * expression it is a type the expression can be evaluated into.
 * @tparam T type to determine plain type of
 */
template <typename T>
struct plain_type<T, require_eigen_t<T>> {
  using type = typename std::decay_t<T>::PlainObject;
};

}  // namespace stan

#endif  // STAN_MATH_PRIM_MAT_META_PLAIN_TYPE_HPP
