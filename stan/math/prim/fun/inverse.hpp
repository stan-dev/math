#ifndef STAN_MATH_PRIM_FUN_INVERSE_HPP
#define STAN_MATH_PRIM_FUN_INVERSE_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the inverse of the specified matrix.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param m specified matrix
 * @return Inverse of the matrix (an empty matrix if the specified matrix has
 * size zero).
 * @throw std::invalid_argument if the matrix is not square.
 */
template <typename EigMat,
          require_eigen_vt<std::is_arithmetic, EigMat>* = nullptr>
inline Eigen::Matrix<value_type_t<EigMat>, EigMat::RowsAtCompileTime,
                     EigMat::ColsAtCompileTime>
inverse(const EigMat& m) {
  check_square("inverse", "m", m);
  if (m.size() == 0) {
    return {};
  }
  return m.inverse();
}

}  // namespace math
}  // namespace stan

#endif
