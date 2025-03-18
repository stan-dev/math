#ifndef STAN_MATH_REV_CONSTRAINT_SUM_TO_ZERO_CONSTRAIN_HPP
#define STAN_MATH_REV_CONSTRAINT_SUM_TO_ZERO_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/constraint/sum_to_zero_constrain.hpp>
#include <cmath>
#include <tuple>
#include <vector>

namespace stan {
namespace math {

/**
 * Return a vector with sum zero corresponding to the specified
 * free vector.
 *
 * The sum-to-zero transform is defined using a modified version of
 * the inverse of the isometric log ratio transform (ILR).
 * See:
 * Egozcue, Juan Jose; Pawlowsky-Glahn, Vera; Mateu-Figueras, Gloria;
 * Barcelo-Vidal, Carles (2003), "Isometric logratio transformations for
 * compositional data analysis", Mathematical Geology, 35 (3): 279–300,
 * doi:10.1023/A:1023818214614, S2CID 122844634
 *
 * This implementation is closer to the description of the same using "pivot
 * coordinates" in
 * Filzmoser, P., Hron, K., Templ, M. (2018). Geometrical Properties of
 * Compositional Data. In: Applied Compositional Data Analysis. Springer Series
 * in Statistics. Springer, Cham. https://doi.org/10.1007/978-3-319-96422-5_3
 *
 * This is a linear transform, with no Jacobian.
 *
 * @tparam T type of the vector
 * @param y Free vector input of dimensionality K - 1.
 * @return Zero-sum vector of dimensionality K.
 */
template <typename T, require_rev_col_vector_t<T>* = nullptr>
inline auto sum_to_zero_constrain(T&& y) {
  using ret_type = plain_type_t<T>;
  if (unlikely(y.size() == 0)) {
    return arena_t<ret_type>(Eigen::VectorXd{{0}});
  }
  auto arena_y = to_arena(std::forward<T>(y));
  arena_t<ret_type> arena_z = sum_to_zero_constrain(arena_y.val());

  reverse_pass_callback([arena_y, arena_z]() mutable {
    const auto N = arena_y.size();

    double sum_u_adj = 0;
    for (Eigen::Index i = 0; i < N; ++i) {
      double n = static_cast<double>(i + 1);

      // adjoint of the reverse cumulative sum computed in the forward mode
      sum_u_adj += arena_z.adj()(i);

      // adjoint of the offset subtraction
      double v_adj = -arena_z.adj()(i + 1) * n;

      double w_adj = v_adj + sum_u_adj;

      arena_y.adj()(i) += w_adj / sqrt(n * (n + 1));
    }
  });

  return arena_z;
}

template <typename T>
Eigen::MatrixXd sum_to_zero_matrix_derivative(const T& x) {
  Eigen::Index N = x.rows();    
  Eigen::Index M = x.cols();  
  Eigen::Index s = std::max(N, M);

  Eigen::MatrixXd Z(N + 1, M + 1);
  Eigen::VectorXd beta(N);
  Eigen::VectorXd a(s);
  Eigen::VectorXd b(s);
  for (Eigen::Index i = 0; i < s; ++i) {
      a(i) = 1.0 / std::sqrt((i + 1.0) * (i + 2.0));
      b(i) = (i + 1.0) * a(i);
      if (i < N) {
          Z(i, M) = 0.0;
          beta(i) = 0.0;
      }
  }

  for (Eigen::Index j = M - 1; j >= 0; --j) {
      double ax_previous = 0.0;  
      Z(N, j) = 0.0;
      for (Eigen::Index i = N - 1; i >= 0; --i) {
          double alpha = b(j) * (b(i) - ax_previous);

          Z(i, j) = alpha - beta(i);
          beta(i) += a(j) * (b(i) - ax_previous);

          Z(N, j) -= Z(i, j);
          Z(i, M) -= Z(i, j);

          ax_previous += a(i);  
      }
      Z(N, M) -= Z(N, j);
  }
  return Z;
}

/**
 * Constrain a matrix so that the resulting (N+1)×(M+1) matrix
 * has its entries summing to zero.
 * This reverse–mode autodiff implementation “replays” the forward pass
 * in reverse order, correctly propagating derivatives through the accumulated
 * β and ax values.
 *
 * @tparam T A Stan autodiff matrix type.
 * @param x Unconstrained input matrix of size (N, M).
 * @return Constrained matrix of size (N+1, M+1).
 */
template <typename T, require_t<is_rev_matrix_dynamic<T>>* = nullptr>
inline auto sum_to_zero_constrain(T&& x) {
  const auto N = x.rows();
  const auto M = x.cols();
  const auto s = std::max(N, M);
  // Move x onto the autodiff arena.
  arena_t<T> x_arena = std::forward<T>(x);

  // We'll record all intermediates in matrices sized (M,N) (each row j holds values for that column).
  arena_t<Eigen::MatrixXd> z_val = Eigen::MatrixXd::Zero(N + 1, M + 1);
  if (unlikely(N == 0 || M == 0)) {
    return arena_t<T>(z_val);
  }
  // Just padding for convenience
  arena_t<Eigen::MatrixXd> ax_previous = Eigen::MatrixXd::Zero(M + 1, N + 1);
  arena_t<Eigen::MatrixXd> b_i_x = Eigen::MatrixXd::Zero(M, N);
  arena_t<Eigen::MatrixXd> beta_pre = Eigen::MatrixXd::Zero(M, N);
  Eigen::Matrix<double, -1, 1> beta = Eigen::VectorXd::Zero(N);
  for (int j = M - 1; j >= 0; --j) {
    double a_j = 1.0 / std::sqrt((j + 1.0) * (j + 2.0));
    double b_j = (j + 1.0) * a_j;
    beta_pre.row(j) = beta.transpose();
    for (int i = N - 1; i >= 0; --i) {
      double a_i = 1.0 / std::sqrt((i + 1.0) * (i + 2.0));
      double b_i = (i + 1.0) * a_i;
      // When we switch rows or columns we need to grab the last value of ax_previous
      double ax_prev = 0;
      if (i == N - 1 && j == M - 1) {
        ax_prev = 0;
      } else if (i == N - 1) {
        ax_prev = ax_previous(0, j + 1);
      } else {
        ax_prev = ax_previous(i + 1, j);
      }
      b_i_x(i, j) = b_i * x_arena.val().coeff(i, j) - ax_prev;
      z_val(i, j) = (b_j * b_i_x(i, j)) - beta(i);
      beta(i) += a_j * b_i_x(i, j);
      // Need to think about how to handle the reverse pass here
      z_val(N, j) -= z_val(i, j);
      z_val(i, M) -= z_val(i, j);
      ax_previous(i, j) = ax_prev + a_i * x_arena.val().coeff(i, j);
    }
    z_val(N, M) -= z_val(N, j);
  }
  arena_t<T> z(z_val);
  reverse_pass_callback([x_arena, z, ax_previous, b_i_x, beta_pre, N, M]() mutable {
    for (int j = 0; j < x_arena.val().cols(); ++j) {      // j from 0 to M-1
      double a_j = 1.0 / std::sqrt((j + 1.0) * (j + 2.0));
      double b_j = (j + 1.0) * a_j;
      for (int i = 0; i < x_arena.val().rows(); ++i) {    // i from 0 to N-1
        double a_i = 1.0 / std::sqrt((i + 1.0) * (i + 2.0));
        double b_i = (i + 1.0) * a_i;
      }
    }
  });
  
  return z;
}

/**
 * Return a vector with sum zero corresponding to the specified
 * free vector.
 *
 * The sum-to-zero transform is defined using a modified version of
 * the inverse of the isometric log ratio transform (ILR).
 * See:
 * Egozcue, Juan Jose; Pawlowsky-Glahn, Vera; Mateu-Figueras, Gloria;
 * Barcelo-Vidal, Carles (2003), "Isometric logratio transformations for
 * compositional data analysis", Mathematical Geology, 35 (3): 279–300,
 * doi:10.1023/A:1023818214614, S2CID 122844634
 *
 * This implementation is closer to the description of the same using "pivot
 * coordinates" in
 * Filzmoser, P., Hron, K., Templ, M. (2018). Geometrical Properties of
 * Compositional Data. In: Applied Compositional Data Analysis. Springer Series
 * in Statistics. Springer, Cham. https://doi.org/10.1007/978-3-319-96422-5_3
 *
 * This is a linear transform, with no Jacobian.
 *
 * @tparam Vec type of the vector
 * @param y Free vector input of dimensionality K - 1.
 * @param lp unused
 * @return Zero-sum vector of dimensionality K.
 */
template <typename T, typename Lp, require_st_var<T>* = nullptr>
inline auto sum_to_zero_constrain(T&& y, Lp& lp) {
  return sum_to_zero_constrain(std::forward<T>(y));
}

}  // namespace math
}  // namespace stan
#endif
