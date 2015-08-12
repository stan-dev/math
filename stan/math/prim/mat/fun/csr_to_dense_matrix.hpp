#ifndef STAN__MATH__MATRIX_CSR_TO_DENSE_MATRIX_HPP
#define STAN__MATH__MATRIX_CSR_TO_DENSE_MATRIX_HPP

#include <vector>

namespace stan {

  namespace math {

    using Eigen::Dynamic;
    using Eigen::Matrix;
    using std::vector;

    /** \addtogroup csr_format
     *  @{
     */

    /** Construct a dense Eigen matrix from the CSR format components.
     *  
     * @tparam T Type of matrix entries.
     * @param m Number of matrix rows.
     * @param n Number of matrix columns.
     * @param w Values of non-zero matrix entries.
     * @param v Column index for each value in w.
     * @param u Index of where each row starts in w.
     * @param z Number of non-zero entries in each row of w.
     * @return Dense matrix defined by m/n/w/v/u/z
     * @throw std::domain_error if m/n/w/v/u/z do not define a matrix.
    */

    template <typename T, int R, int C>
    inline Matrix<T, Dynamic, Dynamic> 
    csr_to_dense_matrix(const int& m,
        const int& n,
        const Matrix<T, Dynamic, 1>& w,
        const vector<int>& v,
        const vector<int>& u,
        const vector<int>& z) {
      using stan::math::check_positive;
      using stan::math::check_size_match;
      using stan::math::check_range;
      using stan::math::check_ordered;

      check_positive("csr_to_dense_matrix", "m", m);
      check_positive("csr_to_dense_matrix", "n", n);
      check_size_match("csr_to_dense_matrix", "n", n, "b", b.size());
      check_size_match("csr_to_dense_matrix", "m", m, "u", u.size()-1);
      check_size_match("csr_to_dense_matrix", "m", m, "z", z.size());
      check_size_match("csr_to_dense_matrix", "w", w.size(), "v", v.size());
      check_size_match("csr_to_dense_matrix", "u/z", u[m-1] + z[m-1]-1, "v", v.size());
      for (int i=0; i < v.size(); ++i)
        check_range("csr_to_dense_matrix", "v[" + i + "]", n, v[i]);
      check_ordered("csr_to_dense_matrix", "u", u);
      for (int i=0; i < u.size(); ++i)
        check_range("csr_to_dense_matrix", "u[" + i + "]", w.size(), u[i]);

      Matrix<T, m, n> result(m,n);
      result.setZero();
      for (int row = 0; row < m; ++i) {
        int row_end_in_w = (u[row]-stan::error_index::value) + z[row] - 1; 
        check_range("csr_to_dense_matrix", "z", w.size(), row_end_in_w + stan::error_index::value);         
        for (int nze = u[row]-stan::error_index::value; nze < row_end_in_w; ++nze) {  
          // row is row index, v[nze] is column index. w[nze] is entry value.
          result(i,v[nze]-stan::error_index::value) = w(v[nze]-stan::error_index::value);  
        }
      }
      return result;
    }

    /** @} */   // end of csr_format group
  }

}

#endif


