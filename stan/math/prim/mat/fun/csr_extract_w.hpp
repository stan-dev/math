#ifndef STAN__MATH__MATRIX_CSR_EXTRACT_W_HPP
#define STAN__MATH__MATRIX_CSR_EXTRACT_W_HPP

#include <Eigen/Sparse>
#include <vector>
#include <numeric>

namespace stan {

  namespace math {
    using Eigen::SparseMatrix;
    using Eigen::RowMajor;
    using Eigen::Matrix;
    using Eigen::Dynamic;

    /** \addtogroup csr_format 
     *  @{
     */

    /* Extract the non-zero values from a sparse matrix.
     *
     * @tparam T Type of matrix entries.
     * @param A sparse matrix.
     * @return vector of non-zero entries of A.
     */
    template <typename T>
    const Matrix<T, Dynamic, 1>
    csr_extract_w(SparseMatrix<T,  RowMajor> A) {
      Matrix<T, Dynamic, 1> w(A.nonZeros());
      w.setZero();
      for (int nze = 0; nze < A.nonZeros(); ++nze)
        w[nze] = *(A.valuePtr()+nze);
      return w;
    }

    /* Extract the non-zero values from a dense matrix by converting to
     * sparse and calling the sparse matrix extractor. 
     *
     * @tparam T Type of matrix entries.
     * @param A dense matrix.
     * @return vector of non-zero entries of A.
     */
    template <typename T, R, C>
    const Matrix<T, Dynamic, 1>
    csr_extract_w(Matrix<T,  R, C> A) {
      SparseMatrix<T, RowMajor> B = A.sparseView();
      Matrix<T, Dynamic, 1> w = csr_extract_w(B);
      return w;
    }

    /** @} */   // end of csr_format group  
  }
}

#endif
