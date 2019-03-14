#ifndef STAN_MATH_PRIM_META_VALUE_TYPE_HPP
#define STAN_MATH_PRIM_META_VALUE_TYPE_HPP

#include <vector>

#include <Eigen/Core>

namespace stan {
namespace math {

/**
 * Primary template class for metaprogram to compute the type of
 * values stored in a container.
 *
 * Only the specializations have behavior that can be used, and
 * all implement a typedef <code>type</code> for the type of the
 * values in the container.
 *
 * tparam T type of container.
 */
template <typename T>
struct value_type {};

/**
 * Template class for metaprogram to compute the type of values
 * stored in a constant container.
 *
 * @tparam T type of container without const modifier.
 */
template <typename T>
struct value_type<const T> {
  typedef typename value_type<T>::type type;
};

/**
 * Template metaprogram class to compute the type of values stored
 * in a standard vector.
 *
 * @tparam T type of elements in standard vector.
 */
template <typename T>
struct value_type<std::vector<T> > {
  /**
   * Type of value stored in a standard vector with type
   * <code>T</code> entries.
   */
  typedef T type;
};

/**
 * Template metaprogram defining the type of values stored in an
 * Eigen matrix, vector, or row vector.
 *
 * @tparam T type of matrix.
 * @tparam R number of rows for matrix.
 * @tparam C number of columns for matrix.
 */
template <typename T, int R, int C>
struct value_type<Eigen::Matrix<T, R, C> > {
  typedef T type;
};

}  // namespace math

}  // namespace stan

#endif
