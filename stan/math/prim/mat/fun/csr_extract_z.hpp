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
    using std::adjacent_difference;

    /** \addtogroup csr_format 
     *  @{
     */

    /* Extract the number of non-zero entries in each row of a sparse
     * matrix.
     *
     * @tparam T Type of matrix entries.
     * @param A sparse matrix.
     * @return vector of counts of non-zero entries in each row of A.
     */
    template <typename T>
    const vector<int>
    csr_extract_z(SparseMatrix<T,  RowMajor> A) {
      vector<int> u(A.outerSize()+1);
      vector<int> z(A.outerSize()+1);
      u = extract_u(A);
      adjacent_difference(u.begin(),  u.end(),  z.begin());
      z.erase(z.begin());
      return z;
    }

    /* Extract the number of non-zero entries in each row of a sparse
     * matrix.
     *
     * @tparam T Type of matrix entries.
     * @param A dense matrix.
     * @return vector of counts of non-zero entries in each row of A.
     */
    template <typename T, R, C>
    const vector<int>
    csr_extract_z(Matrix<T, R, C> A) {
      SparseMatrix<T, RowMajor> B = A.sparseView();
      vector<int> z = csr_extract_z(B);
      return z;
    }
    /** @} */   // end of csr_format group  
  }
}

#endif
