#ifndef STAN__MATH__MATRIX_CSR_EXTRACT_U_HPP
#define STAN__MATH__MATRIX_CSR_EXTRACT_U_HPP

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

    /* Extract the NZE index for each entry from a sparse matrix.
     *
     * @tparam T Type of matrix entries.
     * @param A Sparse matrix.
     * @return vector of indexes into non-zero entries of A.
     */
    template <typename T>
    const vector<int>
    csr_extract_u(SparseMatrix<T,  RowMajor> A) {
      std::vector<int> u(A.outerSize()+1);
      for (int nze = 0; nze <= A.outerSize(); ++nze)
        u[nze] = *(A.outerIndexPtr()+nze) + stan::error_index::value;
      return u;
    }

    /* Extract the NZE index for each entry from a sparse matrix.
     *
     * @tparam T Type of matrix entries.
     * @param A Dense matrix.
     * @return vector of indexes into non-zero entries of A.
     */
    template <typename T, R, C>
    const vector<int>
    csr_extract_u(Matrix<T,  R, C> A) {
      SparseMatrix<T, RowMajor> B = A.sparseView();
      vector<int> u = csr_extract_u(B);
      return u;
    }

    /** @} */   // end of csr_format group  
  }

}

#endif
