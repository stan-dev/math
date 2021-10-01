#ifndef STAN_MATH_PRIM_FUN_CSR_MATRIX_TIMES_VECTOR_HPP
#define STAN_MATH_PRIM_FUN_CSR_MATRIX_TIMES_VECTOR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/csr_u_to_z.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/dot_product.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * @defgroup csr_format Compressed Sparse Row matrix format.
 *  A compressed Sparse Row (CSR) sparse matrix is defined by four
 *  component vectors labeled w, v, and u.  They are defined as:
 *    - w: the non-zero values in the sparse matrix.
 *    - v: column index for each value in w,  as a result this
 *      is the same length as w.
 *    - u: index of where each row starts in w,  length
 *      is equal to the number of rows plus one.  Last entry is
 *      one-past-the-end in w, following the Eigen spec.
 *  Indexing is either zero-based or one-based depending on the value of
 *  stan::error_index::value.  Following the definition of the format in
 *  Eigen, we allow for unused garbage values in w/v which are never read.
 *  All indexing _internal_ to a given function is zero-based.
 *
 *  With only m/n/w/v/u in hand, it is possible to check all
 *  dimensions are sane _except_ the column dimension since it is
 *  implicit.  The error-checking strategy is to check all dimensions
 *  except the column dimension before any work is done inside a
 *  function.  The column index is checked as it is constructed and
 *  used for each entry.  If the column index is not needed it is
 *  not checked.  As a result indexing mistakes might produce non-sensical
 *  operations but out-of-bounds indexing will be caught.
 *
 *  Except for possible garbage values in w/v/u, memory usage is
 *  calculated from the number of non-zero entries (NNZE) and the number of
 *  rows (NR):  2*NNZE + 2*NR + 1.
 */

/**
 * \addtogroup csr_format
 * Return the multiplication of the sparse matrix (specified by
 * by values and indexing) by the specified dense vector.
 *
 * The sparse matrix X of dimension m by n is represented by the
 * vector w (of values), the integer array v (containing one-based
 * column index of each value), the integer array u (containing
 * one-based indexes of where each row starts in w).
 *
 * @tparam T1 type of the sparse matrix
 * @tparam T2 type of the dense vector
 * @param m Number of rows in matrix.
 * @param n Number of columns in matrix.
 * @param w Vector of non-zero values in matrix.
 * @param v Column index of each non-zero value, same
 *          length as w.
 * @param u Index of where each row starts in w, length equal to
 *          the number of rows plus one.
 * @param b Eigen vector which the matrix is multiplied by.
 * @return Dense vector for the product.
 * @throw std::domain_error if m and n are not positive or are nan.
 * @throw std::domain_error if the implied sparse matrix and b are
 *                          not multiplicable.
 * @throw std::invalid_argument if m/n/w/v/u are not internally
 *   consistent, as defined by the indexing scheme.  Extractors are
 *   defined in Stan which guarantee a consistent set of m/n/w/v/u
 *   for a given sparse matrix.
 * @throw std::out_of_range if any of the indexes are out of range.
 */
template <typename T1, typename T2,
          require_all_not_rev_matrix_t<T1, T2>* = nullptr>
inline Eigen::Matrix<return_type_t<T1, T2>, Eigen::Dynamic, 1>
csr_matrix_times_vector(int m, int n, const T1& w, const std::vector<int>& v,
                        const std::vector<int>& u, const T2& b) {
  check_positive("csr_matrix_times_vector", "m", m);
  check_positive("csr_matrix_times_vector", "n", n);
  check_size_match("csr_matrix_times_vector", "n", n, "b", b.size());
  check_size_match("csr_matrix_times_vector", "w", w.size(), "v", v.size());
  check_size_match("csr_matrix_times_vector", "m", m, "u", u.size() - 1);
  check_size_match("csr_matrix_times_vector", "u/z",
                   u[m - 1] + csr_u_to_z(u, m - 1) - 1, "v", v.size());
  for (int i : v) {
    check_range("csr_matrix_times_vector", "v[]", n, i);
  }
  std::vector<int> v_zero_based(v.size());
  std::transform(v.begin(), v.end(), v_zero_based.begin(),
                 [](auto&& x) { return x - 1; });
  std::vector<int> u_zero_based(u.size());
  std::transform(u.begin(), u.end(), u_zero_based.begin(),
                 [](auto&& x) { return x - 1; });
  // u_zero_based[u.size()] = w.size();
  auto&& w_ref = to_ref(w);
  Eigen::Map<const Eigen::SparseMatrix<scalar_type_t<T1>, Eigen::RowMajor>>
      w_mat(m, n, w_ref.size(), u_zero_based.data(), v_zero_based.data(),
            w_ref.data());
  return w_mat * b;
}

}  // namespace math
}  // namespace stan

#endif
