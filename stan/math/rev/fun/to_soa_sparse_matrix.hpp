#ifndef STAN_MATH_REV_FUN_TO_SOA_SPARSE_MATRIX_HPP
#define STAN_MATH_REV_FUN_TO_SOA_SPARSE_MATRIX_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <type_traits>

namespace stan {
namespace math {

/**
 * Create a sparse matrix from the given SoA matrix and indexes.
 * @tparam Options Eigen matrix options.
 * @tparam VarMatrix A @ref var_value with a dense vector inner type
 * @tparam Vec1 Container type of the column indexes.
 * @tparam Vec2 Container type of the row indexes.
 * @param m Number of rows in matrix.
 * @param n Number of columns in matrix.
 * @param w Vector of non-zero values in matrix.
 * @param u Index of where each row starts in w, length equal to
 *          the number of rows plus one.
 * @param v Column index of each non-zero value, same
 *          length as w.
 * @return Sparse matrix.
 */
template <int Options = Eigen::ColMajor, typename VarMatrix, typename Vec1,
          typename Vec2, require_var_t<VarMatrix>* = nullptr,
          require_eigen_dense_base_t<value_type_t<VarMatrix>>* = nullptr,
          require_all_std_vector_vt<std::is_integral, Vec1, Vec2>* = nullptr>
inline auto to_soa_sparse_matrix(int m, int n, VarMatrix&& w, Vec1&& u,
                                 Vec2&& v) {
  auto u_arena = to_arena(std::forward<Vec1>(u));
  auto v_arena = to_arena(std::forward<Vec2>(v));
  using sparse_mat_t = Eigen::SparseMatrix<double, Options>;
  using sparse_arena_mat_t = arena_t<sparse_mat_t>;
  sparse_arena_mat_t arena_val_x(m, n, w.val().size(), u_arena.data(),
                                 v_arena.data(), w.vi_->val_.data());
  sparse_arena_mat_t arena_adj_x(m, n, w.adj().size(), u_arena.data(),
                                 v_arena.data(), w.vi_->adj_.data());
  var_value<sparse_mat_t> var_x(arena_val_x, arena_adj_x);
  return var_x;
}

/**
 * Create a sparse matrix from the given AoS matrix of vars and indexes.
 * @tparam Options Eigen matrix options.
 * @tparam MatrixVar A type inheriting from `Eigen::DenseBase` with a scalar
 * type of @ref var_value
 * @tparam Vec1 Container type of the column indexes.
 * @tparam Vec2 Container type of the row indexes.
 * @param m Number of rows in matrix.
 * @param n Number of columns in matrix.
 * @param w Vector of non-zero values in matrix.
 * @param u Index of where each row starts in w, length equal to
 *          the number of rows plus one.
 * @param v Column index of each non-zero value, same
 *          length as w.
 * @return Sparse matrix.
 */
template <int Options = Eigen::ColMajor, typename MatrixVar, typename Vec1,
          typename Vec2,
          require_eigen_dense_base_vt<is_var, MatrixVar>* = nullptr,
          require_all_std_vector_vt<std::is_integral, Vec1, Vec2>* = nullptr>
inline auto to_soa_sparse_matrix(int m, int n, MatrixVar&& w, Vec1&& u,
                                 Vec2&& v) {
  auto w_arena = to_arena(std::forward<MatrixVar>(w));
  auto u_arena = to_arena(std::forward<Vec1>(u));
  auto v_arena = to_arena(std::forward<Vec2>(v));
  arena_t<Eigen::SparseMatrix<var, Options>> arena_x(
      m, n, w_arena.size(), u_arena.data(), v_arena.data(), w_arena.data());
  var_value<Eigen::SparseMatrix<double, Options>> var_x(value_of(arena_x));
  // No need to copy adj, but need to backprop
  reverse_pass_callback([arena_x, var_x]() mutable {
    using var_sparse_iterator_t =
        typename arena_t<Eigen::SparseMatrix<var, Options>>::InnerIterator;
    using dbl_sparse_iterator_t =
        typename arena_t<Eigen::SparseMatrix<double, Options>>::InnerIterator;
    // arena_x.adj() += var_x.adj() once custom adj() for var sparse matrix
    for (int k = 0; k < arena_x.outerSize(); ++k) {
      var_sparse_iterator_t it_arena_x(arena_x, k);
      dbl_sparse_iterator_t it_var_x(var_x.adj(), k);
      for (; static_cast<bool>(it_arena_x) && static_cast<bool>(it_var_x);
           ++it_arena_x, ++it_var_x) {
        it_arena_x.valueRef().adj() += it_var_x.valueRef();
      }
    }
  });
  return var_x;
}

/**
 * Create a sparse matrix from the given matrix of floats and indexes.
 * @tparam Options Eigen matrix options.
 * @tparam Mat A type inheriting from `Eigen::DenseBase` with an arithmetic
 * scalar type
 * @tparam Vec1 Container type of the column indexes.
 * @tparam Vec2 Container type of the row indexes.
 * @param m Number of rows in matrix.
 * @param n Number of columns in matrix.
 * @param w Vector of non-zero values in matrix.
 * @param u Index of where each row starts in w, length equal to
 *          the number of rows plus one.
 * @param v Column index of each non-zero value, same
 *          length as w.
 * @return Sparse matrix.
 */
template <int Options = Eigen::ColMajor, typename Mat, typename Vec1,
          typename Vec2,
          require_eigen_dense_base_vt<std::is_arithmetic, Mat>* = nullptr,
          require_all_std_vector_vt<std::is_integral, Vec1, Vec2>* = nullptr>
inline auto to_soa_sparse_matrix(int m, int n, Mat&& w, Vec1&& u, Vec2&& v) {
  auto w_arena = to_arena(std::forward<Mat>(w));
  auto u_arena = to_arena(std::forward<Vec1>(u));
  auto v_arena = to_arena(std::forward<Vec2>(v));
  arena_t<Eigen::SparseMatrix<double, Options>> arena_x(
      m, n, w_arena.size(), u_arena.data(), v_arena.data(), w_arena.data());
  return var_value<Eigen::SparseMatrix<double, Options>>(arena_x);
}

}  // namespace math
}  // namespace stan

#endif
