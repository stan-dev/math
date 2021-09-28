#ifndef STAN_MATH_PRIM_FUN_CSR_EXTRACT_U_HPP
#define STAN_MATH_PRIM_FUN_CSR_EXTRACT_U_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <vector>
#include <numeric>

namespace stan {
namespace math {

/** \addtogroup csr_format
 *  @{
 */

/**
 * Extract the NZE index for each entry from a sparse matrix.
 *
 * @tparam T type of elements in the matrix
 * @param A Sparse matrix.
 * @return Array of indexes into non-zero entries of A.
 */
template <typename T>
const std::vector<int> csr_extract_u(
    const Eigen::SparseMatrix<T, Eigen::RowMajor>& A) {
  std::vector<int> u(A.outerSize() + 1);  // last entry is garbage.
  for (int nze = 0; nze <= A.outerSize(); ++nze) {
    u[nze] = *(A.outerIndexPtr() + nze) + stan::error_index::value;
  }
  return u;
}

/**
 * Extract the NZE index for each entry from a sparse matrix.
 *
 * @tparam T type of elements in the matrix
 * @tparam R number of rows, can be Eigen::Dynamic
 * @tparam C number of columns, can be Eigen::Dynamic
 * @param A Dense matrix.
 * @return Array of indexes into non-zero entries of A.
 */
template <typename T, require_eigen_dense_base_t<T>* = nullptr>
const std::vector<int> csr_extract_u(const T& A) {
  // conversion to sparse seems to touch data twice, so we need to call to_ref
  Eigen::SparseMatrix<scalar_type_t<T>, Eigen::RowMajor> B
      = to_ref(A).sparseView();
  std::vector<int> u = csr_extract_u(B);
  return u;
}

/** @} */  // end of csr_format group

}  // namespace math
}  // namespace stan

#endif
