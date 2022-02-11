#ifndef STAN_MATH_PRIM_ERR_IS_MATCHING_DIMS_HPP
#define STAN_MATH_PRIM_ERR_IS_MATCHING_DIMS_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err/is_size_match.hpp>

namespace stan {
namespace math {

/**
 * Return <code>true</code> if the two matrices are of the same size.
 * This function checks the runtime sizes only.
 * @tparam EigMat1 A type derived from `EigenBase`
 * @tparam EigMat2 A type derived from `EigenBase`
 * @param y1 first matrix to test
 * @param y2 second matrix to test
 * @return <code>true</code> if the dimensions of the matrices match
 */
template <typename EigMat1, typename EigMat2,
          require_all_matrix_t<EigMat1, EigMat2>* = nullptr>
inline bool is_matching_dims(const EigMat1& y1, const EigMat2& y2) {
  return is_size_match(y1.rows(), y2.rows())
         && is_size_match(y1.cols(), y2.cols());
}

/**
 * Return <code>true</code> if the two matrices are of the same size.
 * This function checks the runtime sizes and can also check the static
 * sizes as well. For example, a 4x1 matrix is not the same as a vector
 * with 4 elements.
 * @tparam check_compile Whether to check the static sizes
 * @tparam EigMat1 A type derived from `EigenBase`
 * @tparam EigMat2 A type derived from `EigenBase`
 * @param y1 first matrix to test
 * @param y2 second matrix to test
 * @return <code>true</code> if the dimensions of the matrices match
 */
template <bool check_compile, typename EigMat1, typename EigMat2,
          require_all_matrix_t<EigMat1, EigMat2>* = nullptr>
inline bool is_matching_dims(const EigMat1& y1, const EigMat2& y2) {
  return !(check_compile
           && (EigMat1::RowsAtCompileTime != EigMat2::RowsAtCompileTime
               || EigMat1::ColsAtCompileTime != EigMat2::ColsAtCompileTime))
         && is_matching_dims(y1, y2);
}

}  // namespace math
}  // namespace stan
#endif
