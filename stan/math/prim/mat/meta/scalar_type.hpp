#ifndef STAN_MATH_PRIM_MAT_META_SCALAR_TYPE_HPP
#define STAN_MATH_PRIM_MAT_META_SCALAR_TYPE_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/meta/enable_if_eigen.hpp>
#include <stan/math/prim/arr/meta/scalar_type.hpp>

namespace stan {
/**
 * Template metaprogram defining the base scalar type of
 * values stored in an Eigen matrix.
 *
 * @tparam T type of matrix.
 */
template <typename T>
struct scalar_type<T, enable_if_eigen<T>> {
  typedef typename scalar_type<typename std::decay_t<T>::Scalar>::type type;
};

}  // namespace stan
#endif
