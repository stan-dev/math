#ifndef STAN_MATH_PRIM_META_INDEX_TYPE_HPP
#define STAN_MATH_PRIM_META_INDEX_TYPE_HPP

#include <stanh/prim/meta/index_type.hpp>
#include <Eigen/Core>

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
template <typename T, int R, int C>
struct index_type<Eigen::Matrix<T, R, C> > {
  typedef typename Eigen::Matrix<T, R, C>::Index type;
};

}  // namespace math

}  // namespace stan

#endif
#ifndef STAN_MATH_PRIM_META_INDEX_TYPE_HPP
#define STAN_MATH_PRIM_META_INDEX_TYPE_HPP

#include <stanh/prim/meta/index_type.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Template metaprogram class to compute the type of index for a
 * standard vector.
 *
 * @tparam T type of elements in standard vector.
 */
template <typename T>
struct index_type<std::vector<T> > {
  /**
   * Typedef for index of standard vectors.
   */
  typedef typename std::vector<T>::size_type type;
};

}  // namespace math
}  // namespace stan

#endif
#ifndef STAN_MATH_PRIM_META_INDEX_TYPE_HPP
#define STAN_MATH_PRIM_META_INDEX_TYPE_HPP

#include <stanh/prim/meta/index_type.hpp>
#include <Eigen/Core>

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
template <typename T, int R, int C>
struct index_type<Eigen::Matrix<T, R, C> > {
  typedef typename Eigen::Matrix<T, R, C>::Index type;
};

}  // namespace math

}  // namespace stan

#endif
