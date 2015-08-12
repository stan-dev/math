#ifndef STAN__MATH__MATRIX_SPARSE_EXTRACTORS_HPP
#define STAN__MATH__MATRIX_SPARSE_EXTRACTORS_HPP

#include <Eigen/Sparse>
#include <vector>
#include <numeric>

namespace stan {

  namespace math {
    /** \addtogroup csr_format 
     *  @{
     */
    template <typename T>
    const Eigen::Matrix<T, Eigen::Dynamic, 1>
    csr_extract_w(Eigen::SparseMatrix<T,  Eigen::RowMajor> A) {
      Eigen::Matrix<T, Eigen::Dynamic, 1> w(A.nonZeros());
      w.setZero();
      for (int j = 0; j < A.nonZeros(); ++j)
        w[j] = *(A.valuePtr()+j);
      return w;
    }

    template <typename T>
    const std::vector<int>
    csr_extract_v(Eigen::SparseMatrix<T, Eigen::RowMajor> A) {
      std::vector<int> v(A.nonZeros());
      for (int j = 0; j < A.nonZeros(); ++j)
        v[j] = *(A.innerIndexPtr()+j) + 1;  // make 1-indexed
      return v;
    }

    template <typename T>
    const std::vector<int>
    csr_extract_u(Eigen::SparseMatrix<T,  Eigen::RowMajor> A) {
      std::vector<int> u(A.outerSize()+1);
      for (int j = 0; j <= A.outerSize(); ++j)
        u[j] = *(A.outerIndexPtr()+j) + 1;  // make 1-indexed
      return u;
    }

    template <typename T>
    const std::vector<int>
    csr_extract_z(Eigen::SparseMatrix<T,  Eigen::RowMajor> A) {
      std::vector<int> u(A.outerSize()+1);
      std::vector<int> z(A.outerSize()+1);
      u = extract_u(A);
      std::adjacent_difference(u.begin(),  u.end(),  z.begin());
      z.erase(z.begin());
      return z;
    }

    /** @} */   // end of csr_format group  
  }
}

#endif
