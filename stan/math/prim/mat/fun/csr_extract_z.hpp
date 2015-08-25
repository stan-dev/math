#ifndef STAN_MATH_PRIM_MAT_FUN_CSR_EXTRACT_Z_HPP
#define STAN_MATH_PRIM_MAT_FUN_CSR_EXTRACT_Z_HPP

#include <Eigen/Sparse>
#include <vector>
#include <numeric>

namespace stan {

  namespace math {
    using Eigen::SparseMatrix;
    using Eigen::RowMajor;
    using Eigen::Matrix;
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
    csr_extract_z(const SparseMatrix<T,  RowMajor>& A) {
      vector<int> u(A.outerSize()+1);
      vector<int> z(A.outerSize()+1);
      u = csr_extract_u(A);
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
    template <typename T, int R, int C>
    const vector<int>
    csr_extract_z(const Matrix<T, R, C>& A) {
      SparseMatrix<T, RowMajor> B = A.sparseView();
      vector<int> z = csr_extract_z(B);
      return z;
    }
    /** @} */   // end of csr_format group
  }
}

#endif
