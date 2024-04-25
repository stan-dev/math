#ifndef STAN_MATH_REV_FUN_CSR_MATRIX_TIMES_VECTOR_HPP
#define STAN_MATH_REV_FUN_CSR_MATRIX_TIMES_VECTOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/fun/to_soa_sparse_matrix.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/csr_u_to_z.hpp>
#include <vector>

namespace stan {
namespace math {

namespace internal {
/**
 * `vari` for csr_matrix_times_vector
 * @note `csr_matrix_times_vector` uses the old inheritance
 *  style to set up the reverse pass because of a linking
 *  issue on windows when using flto.
 *
 * @tparam Result_ Either a type inheriting from `Eigen::DenseBase` with scalar
 * type `var` or a `var<T>` where `T` inherits from `Eigen::DenseBase`
 * @tparam WMat_ Either a type inheriting from `Eigen::DenseBase` with scalar
 * type `var` or `double`. Or a `var<T>` where `T` inherits from
 * `Eigen::SparseBase`
 * @tparam B_ Either a type inheriting from `Eigen::DenseBase` with scalar type
 * `var` or `double`. Or a `var<T>` where `T` inherits from `Eigen::DenseBase`
 *
 */
template <typename Result_, typename WMat_, typename B_>
struct csr_adjoint : public vari {
  std::decay_t<Result_> res_;
  std::decay_t<WMat_> w_mat_;
  std::decay_t<B_> b_;

  template <typename T1, typename T2, typename T3>
  csr_adjoint(T1&& res, T2&& w_mat, T3&& b)
      : vari(0.0),
        res_(std::forward<T1>(res)),
        w_mat_(std::forward<T2>(w_mat)),
        b_(std::forward<T3>(b)) {}

  void chain() { chain_internal(res_, w_mat_, b_); }

  /**
   * Overload for calculating adjoints of `w_mat` and `b`
   * @tparam Result Either a type inheriting from `Eigen::DenseBase` with scalar
   * type `var` or a `var<T>` where `T` inherits from `Eigen::DenseBase`
   * @tparam WMat Either a type inheriting from `Eigen::DenseBase` with scalar
   * type `var`. Or a `var<T>` where `T` inherits from `Eigen::SparseBase`
   * @tparam B Either a type inheriting from `Eigen::DenseBase` with scalar type
   * `var`. Or a `var<T>` where `T` inherits from `Eigen::DenseBase`
   * @param res The vector result of the forward pass calculation
   * @param w_mat A sparse matrix
   * @param b A vector
   */
  template <typename Result, typename WMat, typename B,
            require_rev_matrix_t<WMat>* = nullptr,
            require_rev_matrix_t<B>* = nullptr>
  inline void chain_internal(Result&& res, WMat&& w_mat, B&& b) {
    w_mat.adj() += res.adj() * b.val().transpose();
    b.adj() += w_mat.val().transpose() * res.adj();
  }

  /**
   * Overload for calculating adjoints of `w_mat`
   * @tparam Result Either a type inheriting from `Eigen::DenseBase` with scalar
   * type `var` or a `var<T>` where `T` inherits from `Eigen::DenseBase`
   * @tparam WMat Either a type inheriting from `Eigen::DenseBase` with scalar
   * type `var`. Or a `var<T>` where `T` inherits from `Eigen::SparseBase`
   * @tparam B Either a type inheriting from `Eigen::DenseBase` with scalar type
   * `double`
   * @param res The vector result of the forward pass calculation
   * @param w_mat A sparse matrix
   * @param b A vector
   */
  template <typename Result, typename WMat, typename B,
            require_rev_matrix_t<WMat>* = nullptr,
            require_not_rev_matrix_t<B>* = nullptr>
  inline void chain_internal(Result&& res, WMat&& w_mat, B&& b) {
    w_mat.adj() += res.adj() * b.transpose();
  }

  /**
   * Overload for calculating adjoints of `b`
   * @tparam Result Either a type inheriting from `Eigen::DenseBase` with scalar
   * type `var` or a `var<T>` where `T` inherits from `Eigen::DenseBase`
   * @tparam WMat Either a type inheriting from `Eigen::DenseBase` with scalar
   * type `double`
   * @tparam B Either a type inheriting from `Eigen::DenseBase` with scalar type
   * `var` or a `var<T>` where `T` inherits from `Eigen::DenseBase`
   * @param res The vector result of the forward pass calculation
   * @param w_mat A sparse matrix
   * @param b A vector
   */
  template <typename Result, typename WMat, typename B,
            require_not_rev_matrix_t<WMat>* = nullptr,
            require_rev_matrix_t<B>* = nullptr>
  inline void chain_internal(Result&& res, WMat&& w_mat, B&& b) {
    b.adj() += w_mat.transpose() * res.adj();
  }
};

/**
 * Helper function to construct the csr_adjoint struct.
 * @tparam Result_ Either a type inheriting from `Eigen::DenseBase` with scalar
 * type `var` or a `var<T>` where `T` inherits from `Eigen::DenseBase`
 * @tparam WMat_ Either a type inheriting from `Eigen::DenseBase` with scalar
 * type `var` or `double`. Or a `var<T>` where `T` inherits from
 * `Eigen::SparseBase`
 * @tparam B_ Either a type inheriting from `Eigen::DenseBase` with scalar type
 * `var` or `double`. Or a `var<T>` where `T` inherits from `Eigen::DenseBase`
 *
 * @param res The vector result of the forward pass calculation
 * @param w_mat A sparse matrix
 * @param b A vector
 */
template <typename Result_, typename WMat_, typename B_>
inline void make_csr_adjoint(Result_&& res, WMat_&& w_mat, B_&& b) {
  new csr_adjoint<std::decay_t<Result_>, std::decay_t<WMat_>, std::decay_t<B_>>(
      std::forward<Result_>(res), std::forward<WMat_>(w_mat),
      std::forward<B_>(b));
  return;
}
}  // namespace internal

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
template <typename T1, typename T2, require_any_rev_matrix_t<T1, T2>* = nullptr>
inline auto csr_matrix_times_vector(int m, int n, const T1& w,
                                    const std::vector<int>& v,
                                    const std::vector<int>& u, const T2& b) {
  using sparse_val_mat
      = Eigen::Map<const Eigen::SparseMatrix<double, Eigen::RowMajor>>;
  using sparse_dense_mul_type
      = decltype((std::declval<sparse_val_mat>() * value_of(b)).eval());
  using return_t = return_var_matrix_t<sparse_dense_mul_type, T1, T2>;

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
  std::vector<int, arena_allocator<int>> v_arena(v.size());
  std::transform(v.begin(), v.end(), v_arena.begin(),
                 [](auto&& x) { return x - 1; });
  std::vector<int, arena_allocator<int>> u_arena(u.size());
  std::transform(u.begin(), u.end(), u_arena.begin(),
                 [](auto&& x) { return x - 1; });
  using sparse_var_value_t
      = var_value<Eigen::SparseMatrix<double, Eigen::RowMajor>>;
  if (!is_constant<T2>::value && !is_constant<T1>::value) {
    arena_t<promote_scalar_t<var, T2>> b_arena = b;
    sparse_var_value_t w_mat_arena
        = to_soa_sparse_matrix<Eigen::RowMajor>(m, n, w, u_arena, v_arena);
    arena_t<return_t> res = w_mat_arena.val() * value_of(b_arena);
    stan::math::internal::make_csr_adjoint(res, w_mat_arena, b_arena);
    return return_t(res);
  } else if (!is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T2>> b_arena = b;
    auto w_val_arena = to_arena(value_of(w));
    sparse_val_mat w_val_mat(m, n, w_val_arena.size(), u_arena.data(),
                             v_arena.data(), w_val_arena.data());
    arena_t<return_t> res = w_val_mat * value_of(b_arena);
    stan::math::internal::make_csr_adjoint(res, w_val_mat, b_arena);
    return return_t(res);
  } else {
    sparse_var_value_t w_mat_arena
        = to_soa_sparse_matrix<Eigen::RowMajor>(m, n, w, u_arena, v_arena);
    auto b_arena = to_arena(value_of(b));
    arena_t<return_t> res = w_mat_arena.val() * b_arena;
    stan::math::internal::make_csr_adjoint(res, w_mat_arena, b_arena);
    return return_t(res);
  }
}

}  // namespace math
}  // namespace stan

#endif
