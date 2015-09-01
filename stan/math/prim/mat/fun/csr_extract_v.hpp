#ifndef STAN_MATH_PRIM_MAT_FUN_CSR_EXTRACT_V_HPP
#define STAN_MATH_PRIM_MAT_FUN_CSR_EXTRACT_V_HPP

#include <stan/math.hpp>
#include <Eigen/Sparse>
#include <vector>
#include <numeric>

namespace stan {

  namespace math {
    using Eigen::SparseMatrix;
    using Eigen::RowMajor;
    using Eigen::Matrix;
    using std::vector;

    /** \addtogroup csr_format 
     *  @{
     */

    /* Extract the column indexes for non-zero value from a sparse
     * matrix.
     * @tparam T Type of matrix entries.
     * @param A sparse matrix.
     * @return vector of column indexes for non-zero entries of A.
     */
    template <typename T>
    const vector<int>
    csr_extract_v(const SparseMatrix<T, RowMajor>& A) {
      vector<int> v(A.nonZeros());
      for (int nze = 0; nze < A.nonZeros(); ++nze)
        v[nze] = *(A.innerIndexPtr()+nze) + stan::error_index::value;
      return v;
    }

    /* Extract the column indexes for non-zero values from a dense 
     * matrix by converting to sparse and calling the sparse matrix 
     * extractor. 
     *
     * @tparam T Type of matrix entries.
     * @param A dense matrix.
     * @return vector of column indexes to non-zero entries of A.
     */
    template <typename T, int R, int C>
    const vector<int>
    csr_extract_v(const Matrix<T,  R, C>& A) {
      SparseMatrix<T, RowMajor> B = A.sparseView();
      vector<int> v = csr_extract_v(B);
      return v;
    }

    /** @} */   // end of csr_format group
  }
}

#endif
