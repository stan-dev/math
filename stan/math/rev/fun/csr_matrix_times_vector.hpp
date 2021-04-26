#ifndef STAN_MATH_REV_FUN_CSR_MATRIX_TIMES_VECTOR_HPP
#define STAN_MATH_REV_FUN_CSR_MATRIX_TIMES_VECTOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/csr_u_to_z.hpp>
#include <vector>

namespace stan {
namespace math {

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
template <typename T1, typename T2, require_st_arithmetic<T1>* = nullptr,
          require_st_var<T2>* = nullptr>
inline auto csr_matrix_times_vector(int m, int n, const T1& w,
                                    const std::vector<int>& v,
                                    const std::vector<int>& u, const T2& b) {
  check_positive("csr_matrix_times_vector", "m", m);
  check_positive("csr_matrix_times_vector", "n", n);
  check_size_match("csr_matrix_times_vector", "n", n, "b", b.size());
  check_size_match("csr_matrix_times_vector", "m", m, "u", u.size() - 1);
  check_size_match("csr_matrix_times_vector", "w", w.size(), "v", v.size());
  check_size_match("csr_matrix_times_vector", "u/z",
                   u[m - 1] + csr_u_to_z(u, m - 1) - 1, "v", v.size());
  for (int i : v) {
    check_range("csr_matrix_times_vector", "v[]", n, i);
  }
  std::vector<int, arena_allocator<int>> arena_v(v.begin(), v.end());
  std::vector<int, arena_allocator<int>> arena_u(u.begin(), u.end());
  auto arena_w = to_arena(w);
  auto arena_b = to_arena(b);
  Eigen::Map<Eigen::SparseMatrix<scalar_type_t<T1>>> arena_sp_map(
      m, n, arena_b.size(), arena_v.data(), arena_u.data(), arena_w.data());
  using sparse_dense_mul_type
      = decltype((arena_sp_map * value_of(arena_b)).eval());
  using return_t = return_var_matrix_t<sparse_dense_mul_type, T1, T2>;
  arena_t<return_t> result = arena_sp_map * arena_b.val();
  reverse_pass_callback([arena_v, arena_u, arena_w, arena_b, result, m,
                         n]() mutable {
    Eigen::Map<Eigen::SparseMatrix<scalar_type_t<T1>>> arena_sp_map(
        m, n, arena_b.size(), arena_v.data(), arena_u.data(), arena_w.data());
    arena_b.adj() += arena_sp_map.transpose() * result.adj();
  });

  return return_t(result);
}

}  // namespace math
}  // namespace stan

#endif
