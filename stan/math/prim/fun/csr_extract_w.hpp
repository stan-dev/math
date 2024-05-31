#ifndef STAN_MATH_PRIM_FUN_CSR_EXTRACT_W_HPP
#define STAN_MATH_PRIM_FUN_CSR_EXTRACT_W_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>

namespace stan {
namespace math {

/** \addtogroup csr_format
 *  @{
 */

/* Extract the non-zero values from a sparse matrix.
 *
 * @tparam T type of elements in the matrix
 * @param[in] A sparse matrix.
 * @return Vector of non-zero entries of A.
 */
template <typename T>
const Eigen::Matrix<T, Eigen::Dynamic, 1> csr_extract_w(
    const Eigen::SparseMatrix<T, Eigen::RowMajor>& A) {
  auto a_nonzeros = A.nonZeros();
  Eigen::Matrix<T, Eigen::Dynamic, 1> w
      = Eigen::Matrix<T, Eigen::Dynamic, 1>::Zero(a_nonzeros);
  for (int nze = 0; nze < a_nonzeros; ++nze) {
    w[nze] = *(A.valuePtr() + nze);
  }
  return w;
}

/* Extract the non-zero values from a dense matrix by converting
 * to sparse and calling the sparse matrix extractor.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 *
 * @param[in] A dense matrix.
 * @return Vector of non-zero entries of A.
 */
template <typename T, require_eigen_dense_base_t<T>* = nullptr>
const Eigen::Matrix<scalar_type_t<T>, Eigen::Dynamic, 1> csr_extract_w(
    const T& A) {
  // conversion to sparse seems to touch data twice, so we need to call to_ref
  Eigen::SparseMatrix<scalar_type_t<T>, Eigen::RowMajor> B
      = to_ref(A).sparseView();
  return csr_extract_w(B);
}

/** @} */  // end of csr_format group

}  // namespace math
}  // namespace stan

#endif
