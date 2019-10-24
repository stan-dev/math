#ifndef STAN_MATH_PRIM_MAT_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_MAT_META_SCALAR_TYPE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/is_eigen.hpp>
#include <stan/math/prim/arr/meta/scalar_type.hpp>
#include <type_traits>
namespace stan {
/**
 * Template metaprogram defining the base scalar type of
 * values stored in an Eigen matrix.
 *
 * @tparam T type of matrix.
 */
template <typename T>
struct scalar_type<T, std::enable_if_t<is_eigen<T>::value>> {
  using type = scalar_type_t<typename std::decay_t<T>::Scalar>;
};
}  // namespace stan
#endif
