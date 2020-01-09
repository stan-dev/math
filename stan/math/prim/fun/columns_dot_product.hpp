#ifndef STAN_MATH_PRIM_FUN_COLUMNS_DOT_PRODUCT_HPP
#define STAN_MATH_PRIM_FUN_COLUMNS_DOT_PRODUCT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of the specified vectors.
 *
 * @tparam R1 number of rows, can be Eigen::Dynamic
 * @tparam C1 number of columns, can be Eigen::Dynamic
 * @tparam R2 number of rows, can be Eigen::Dynamic
 * @tparam C2 number of columns, can be Eigen::Dynamic
 *
 * @param v1 First vector.
 * @param v2 Second vector.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <int R1, int C1, int R2, int C2>
inline Eigen::Matrix<double, 1, C1> columns_dot_product(
    const Eigen::Matrix<double, R1, C1>& v1,
    const Eigen::Matrix<double, R2, C2>& v2) {
  check_matching_sizes("columns_dot_product", "v1", v1, "v2", v2);
  return (v1.transpose() * v2).diagonal();
}

}  // namespace math
}  // namespace stan

#endif
