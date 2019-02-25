#ifndef STAN_MATH_PRIM_META_VALUE_TYPE_HPP
#define STAN_MATH_PRIM_META_VALUE_TYPE_HPP

#include <stanh/prim/meta/value_type.hpp>
#include <Eigen/Core>

namespace stan {
namespace math {

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
#ifndef STAN_MATH_PRIM_META_VALUE_TYPE_HPP
#define STAN_MATH_PRIM_META_VALUE_TYPE_HPP

#include <stanh/prim/meta/value_type.hpp>
#include <vector>

namespace stan {
namespace math {

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

}  // namespace math
}  // namespace stan
#endif
#ifndef STAN_MATH_PRIM_META_VALUE_TYPE_HPP
#define STAN_MATH_PRIM_META_VALUE_TYPE_HPP

#include <stanh/prim/meta/value_type.hpp>
#include <Eigen/Core>

namespace stan {
namespace math {

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
