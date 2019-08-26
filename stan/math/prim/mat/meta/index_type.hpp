#ifndef STAN_MATH_PRIM_MAT_META_INDEX_TYPE_HPP
#define STAN_MATH_PRIM_MAT_META_INDEX_TYPE_HPP

#include <stan/math/prim/scal/meta/index_type.hpp>
#include <stan/math/prim/mat/meta/is_eigen.hpp>
#include <Eigen/Core>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Template metaprogram defining typedef for the type of index for
 * an Eigen matrix, vector, or row vector.
 *
 * @tparam T type of matrix.
 * @tparam R number of rows for matrix.
 * @tparam C number of columns for matrix.
 */
template <typename T>
struct index_type<T, std::enable_if_t<is_eigen<T>::value>> {
  typedef typename std::decay_t<T>::Index type;
};

}  // namespace math

}  // namespace stan

#endif
