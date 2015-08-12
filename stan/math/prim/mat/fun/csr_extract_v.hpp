#ifndef STAN__MATH__MATRIX_SPARSE_EXTRACTORS_HPP
#define STAN__MATH__MATRIX_SPARSE_EXTRACTORS_HPP

#include <Eigen/Sparse>
#include <vector>
#include <numeric>

namespace stan {

  namespace math {
    using Eigen::SparseMatrix;
    using Eigen::RowMajor;
    using Eigen::Matrix;
    using Eigen::Dynamic;
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
    csr_extract_v(SparseMatrix<T, RowMajor> A) {
      vector<int> v(A.nonZeros());
      for (int nze = 0; nze < A.nonZeros(); ++j)
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
    template <typename T, R, C>
    const vector<int> 
    csr_extract_w(Matrix<T,  R, C> A) {
      SparseMatrix<T, RowMajor> B = A.sparseView();
      vector<int> v = csr_extract_v(B);
      return v;
    }

    /** @} */   // end of csr_format group  
  }
}

#endif
