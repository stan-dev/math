#include <stan/math/prim/scal/meta/broadcast_array.hpp>
#include <Eigen/Dense>
#include <stdexcept>

#ifndef STAN_MATH_PRIM_MAT_META_ASSIGN_TO_MATRIX_OR_BROADCAST_ARRAY_HPP
#define STAN_MATH_PRIM_MAT_META_ASSIGN_TO_MATRIX_OR_BROADCAST_ARRAY_HPP

namespace stan {
namespace math {
/**
 * This program is used to assign to either
 * an Eigen::Matrix or to a stan::math::internal::broadcast_array.
 *
 * @tparam T Type of argument to be assigned to.
 * @tparam S Type of argument of what we are assigning.
 * @param arg1 Argument to be assigned to.
 * @param arg2 Argument we are assigning.
 */
template <typename T, typename S>
void assign_to_matrix_or_broadcast_array(T &arg1, S &arg2);

/**
 * This program is used to assign to a
 * stan::math::internal::broadcast_array.
 *
 * @tparam TT Type of elements of broadcast_array we are assigning to.
 * @tparam S Type of argument of what we are assigning.
 * @param arg1 Argument to be assigned to.
 * @param arg2 Argument we are assigning.
 */
template <typename TT, typename S>
void assign_to_matrix_or_broadcast_array(internal::broadcast_array<TT> &arg1,
                                         S &arg2) {
  arg1[0] = arg2(0, 0);
}

/**
 * This program is used to assign to an Eigen::Matrix.
 *
 * @tparam TT Type of entries of Eigen::Matrix we are assigning to.
 * @tparam R Number of rows of this matrix at compile time.
 * @tparam C Number of columns of this matrix at compile time.
 * @tparam S Type of argument of what we are assigning.
 * @param arg1 Argument to be assigned to.
 * @param arg2 Argument we are assigning.
 */
template <typename TT, int R, int C, typename S>
void assign_to_matrix_or_broadcast_array(Eigen::Matrix<TT, R, C> &arg1,
                                         S &arg2) {
  arg1 = arg2;
}
}  // namespace math
}  // namespace stan

#endif
