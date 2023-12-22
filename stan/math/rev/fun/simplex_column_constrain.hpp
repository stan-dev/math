#ifndef STAN_MATH_REV_FUN_SIMPLEX_COLUMN_CONSTRAIN_HPP
#define STAN_MATH_REV_FUN_SIMPLEX_COLUMN_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/fun/value_of.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv_logit.hpp>
#include <stan/math/prim/fun/log1p_exp.hpp>
#include <cmath>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

namespace internal {
template <typename Mat, int NewOptions, typename = void>
struct change_eigen_options_impl {};

template <typename Mat, int NewOptions>
struct change_eigen_options_impl<Mat, NewOptions,
                                 require_eigen_matrix_base_t<Mat>> {
  using type
      = Eigen::Matrix<typename Mat::Scalar, Mat::RowsAtCompileTime,
                      Mat::ColsAtCompileTime, NewOptions,
                      Mat::MaxRowsAtCompileTime, Mat::MaxColsAtCompileTime>;
};

template <typename Mat, int NewOptions>
struct change_eigen_options_impl<var_value<Mat>, NewOptions,
                                 require_eigen_matrix_base_t<Mat>> {
  using type = var_value<Eigen::Matrix<
      typename Mat::Scalar, Mat::RowsAtCompileTime, Mat::ColsAtCompileTime,
      NewOptions, Mat::MaxRowsAtCompileTime, Mat::MaxColsAtCompileTime>>;
};

template <typename Mat, int NewOptions>
struct change_eigen_options_impl<Mat, NewOptions, require_eigen_array_t<Mat>> {
  using type
      = Eigen::Array<typename Mat::Scalar, Mat::RowsAtCompileTime,
                     Mat::ColsAtCompileTime, NewOptions,
                     Mat::MaxRowsAtCompileTime, Mat::MaxColsAtCompileTime>;
};

template <typename Mat, int NewOptions>
using change_eigen_options_t =
    typename change_eigen_options_impl<plain_type_t<std::decay_t<Mat>>,
                                       NewOptions>::type;
}  // namespace internal

/**
 * Return a matrix with columns as simplex vectors.
 * A simplex is a vector containing values greater than or equal
 * to 0 that sum to 1.  A matrix (K-1, M) unconstrained values
 * will produce a matrix of simplex vectors of size (K, M).
 *
 * The transform is based on a centered stick-breaking process.
 *
 * @tparam T Type of matrix to constrain
 * @param y Free matrix input of dimensionality (K - 1, M)
 * @return matrix of column simplexes of dimensionality (K, M)
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto simplex_column_constrain(const T& y) {
  using ret_type = plain_type_t<T>;
  const Eigen::Index N = y.rows();
  const Eigen::Index M = y.cols();
  using eigen_mat_rowmajor
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  arena_t<eigen_mat_rowmajor> x_val(N + 1, M);
  if (unlikely(N == 0 || M == 0)) {
    return ret_type(x_val);
  }
  arena_t<internal::change_eigen_options_t<T, Eigen::RowMajor>> arena_y = y;
  arena_t<eigen_mat_rowmajor> arena_z(N, M);
  using arr_vec = Eigen::Array<double, 1, -1>;
  arr_vec stick_len = arr_vec::Constant(M, 1.0);
  for (Eigen::Index k = 0; k < N; ++k) {
    const double log_N_minus_k = std::log(N - k);
    arena_z.row(k)
        = inv_logit(arena_y.array().row(k).val_op() - log_N_minus_k).matrix();
    x_val.row(k) = stick_len.array() * arena_z.array().row(k);
    stick_len -= x_val.array().row(k);
  }
  x_val.row(N) = stick_len;
  arena_t<ret_type> arena_x = x_val;
  reverse_pass_callback([arena_y, arena_x, arena_z]() mutable {
    const Eigen::Index N = arena_y.rows();
    auto arena_x_arr = arena_x.array();
    auto arena_y_arr = arena_y.array();
    auto arena_z_arr = arena_z.array();
    auto stick_len_val = arena_x.array().row(N).val().eval();
    auto stick_len_adj = arena_x.array().row(N).adj().eval();
    for (Eigen::Index k = N; k-- > 0;) {
      arena_x_arr.row(k).adj() -= stick_len_adj;
      stick_len_val += arena_x_arr.row(k).val();
      stick_len_adj += arena_x_arr.row(k).adj() * arena_z_arr.row(k);
      auto arena_z_adj = arena_x_arr.row(k).adj() * stick_len_val;
      arena_y_arr.row(k).adj()
          += arena_z_adj * arena_z_arr.row(k) * (1.0 - arena_z_arr.row(k));
    }
  });
  return ret_type(arena_x);
}

/**
 * Return a matrix with columns as simplex vectors
 * and increment the specified log probability reference with
 * the log absolute Jacobian determinant of the transform.
 *
 * The simplex transform is defined through a centered
 * stick-breaking process.
 *
 * @tparam T type of the matrix to constrain
 * @param y Free matrix input of dimensionality N, M.
 * @param lp Log probability reference to increment.
 * @return Matrix of simplex columns of dimensionality (N + 1, M).
 */
template <typename T, require_rev_matrix_t<T>* = nullptr>
auto simplex_column_constrain(const T& y, scalar_type_t<T>& lp) {
  using ret_type = plain_type_t<T>;
  const Eigen::Index N = y.rows();
  const Eigen::Index M = y.cols();
  using eigen_mat_rowmajor
      = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  arena_t<eigen_mat_rowmajor> x_val(N + 1, M);
  if (unlikely(N == 0 || M == 0)) {
    return ret_type(x_val);
  }
  arena_t<internal::change_eigen_options_t<T, Eigen::RowMajor>> arena_y = y;
  arena_t<eigen_mat_rowmajor> arena_z(N, M);
  using arr_vec = Eigen::Array<double, 1, -1>;
  arr_vec stick_len = arr_vec::Constant(M, 1.0);
  arr_vec adj_y_k(N);
  for (Eigen::Index k = 0; k < N; ++k) {
    double log_N_minus_k = std::log(N - k);
    adj_y_k = arena_y.array().row(k).val() - log_N_minus_k;
    arena_z.array().row(k) = inv_logit(adj_y_k);
    x_val.array().row(k) = stick_len * arena_z.array().row(k);
    lp += sum(log(stick_len)) - sum(log1p_exp(-adj_y_k))
          - sum(log1p_exp(adj_y_k));
    stick_len -= x_val.array().row(k);
  }
  x_val.array().row(N) = stick_len;
  arena_t<ret_type> arena_x = x_val;
  reverse_pass_callback([arena_y, arena_x, arena_z, lp]() mutable {
    const Eigen::Index N = arena_y.rows();
    auto arena_x_arr = arena_x.array();
    auto arena_y_arr = arena_y.array();
    auto arena_z_arr = arena_z.array();
    auto stick_len_val = arena_x.array().row(N).val().eval();
    auto stick_len_adj = arena_x.array().row(N).adj().eval();
    for (Eigen::Index k = N; k-- > 0;) {
      const double log_N_minus_k = std::log(N - k);
      arena_x_arr.row(k).adj() -= stick_len_adj;
      stick_len_val += arena_x_arr.row(k).val();
      stick_len_adj += lp.adj() / stick_len_val
                       + arena_x_arr.row(k).adj() * arena_z_arr.row(k);
      auto adj_y_k = arena_y_arr.row(k).val() - log_N_minus_k;
      auto arena_z_adj = arena_x_arr.row(k).adj() * stick_len_val;
      arena_y_arr.row(k).adj()
          += -(lp.adj() * inv_logit(adj_y_k)) + lp.adj() * inv_logit(-adj_y_k)
             + arena_z_adj * arena_z_arr.row(k) * (1.0 - arena_z_arr.row(k));
    }
  });
  return ret_type(arena_x);
}

}  // namespace math
}  // namespace stan
#endif
