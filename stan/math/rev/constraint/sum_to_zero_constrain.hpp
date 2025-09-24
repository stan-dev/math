#ifndef STAN_MATH_REV_CONSTRAINT_SUM_TO_ZERO_CONSTRAIN_HPP
#define STAN_MATH_REV_CONSTRAINT_SUM_TO_ZERO_CONSTRAIN_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/inv_sqrt.hpp>
#include <stan/math/prim/fun/sqrt.hpp>
#include <stan/math/prim/constraint/sum_to_zero_constrain.hpp>
#include <cmath>

namespace stan {
namespace math {

namespace internal {

/**
 * The reverse pass backprop for the sum_to_zero_constrain on
 * vectors. This is separated out so it can also be called by
 * simplex_constrain.
 *
 * @tparam T type of the adjoint vector
 * @param y_adj The adjoint of the free vector (size N)
 * @param z_adj The adjoint of the zero-sum vector (size N + 1)
 */
template <typename T>
void sum_to_zero_vector_backprop(T&& y_adj, const Eigen::VectorXd& z_adj) {
  const auto N = y_adj.size();

  double sum_u_adj = 0;
  for (int i = 0; i < N; ++i) {
    double n = static_cast<double>(i + 1);

    // adjoint of the reverse cumulative sum computed in the forward mode
    sum_u_adj += z_adj.coeff(i);

    // adjoint of the offset subtraction
    double v_adj = -z_adj.coeff(i + 1) * n;

    double w_adj = v_adj + sum_u_adj;

    y_adj.coeffRef(i) += w_adj / sqrt(n * (n + 1));
  }
}

}  // namespace internal

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
    internal::sum_to_zero_vector_backprop(arena_y.adj(), arena_z.adj());
  });

  return arena_z;
}

/**
 * Return a matrix that sums to zero over both the rows
 * and columns corresponding to the free matrix x.
 *
 * This is a linear transform, with no Jacobian.
 *
 * @tparam Mat type of the matrix
 * @param x Free matrix input of dimensionality (N - 1, M - 1).
 * @return Zero-sum matrix of dimensionality (N, M).
 */
template <typename T, require_rev_matrix_t<T>* = nullptr,
          require_not_t<is_rev_vector<T>>* = nullptr>
inline auto sum_to_zero_constrain(T&& x) {
  using ret_type = plain_type_t<T>;
  if (unlikely(x.size() == 0)) {
    return arena_t<ret_type>(Eigen::MatrixXd{{0}});
  }
  auto arena_x = to_arena(std::forward<T>(x));
  arena_t<ret_type> arena_z = sum_to_zero_constrain(arena_x.val());

  reverse_pass_callback([arena_x, arena_z]() mutable {
    const auto Nf = arena_x.val().rows();
    const auto Mf = arena_x.val().cols();

    Eigen::VectorXd d_beta = Eigen::VectorXd::Zero(Nf);

    for (int j = 0; j < Mf; ++j) {
      double a_j = inv_sqrt((j + 1.0) * (j + 2.0));
      double b_j = (j + 1.0) * a_j;

      double d_ax = 0.0;

      for (int i = 0; i < Nf; ++i) {
        double a_i = inv_sqrt((i + 1.0) * (i + 2.0));
        double b_i = (i + 1.0) * a_i;

        double dY = arena_z.adj().coeff(i, j) - arena_z.adj().coeff(Nf, j)
                    + arena_z.adj().coeff(Nf, Mf) - arena_z.adj().coeff(i, Mf);
        double dI_from_beta = a_j * d_beta.coeff(i);
        d_beta.coeffRef(i) += -dY;

        double dI_from_alpha = b_j * dY;
        double dI = dI_from_alpha + dI_from_beta;
        arena_x.adj().coeffRef(i, j) += b_i * dI + a_i * d_ax;
        d_ax -= dI;
      }
    }
  });

  return arena_z;
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
 * @tparam T type of the vector or matrix
 * @param y Free vector or matrix.
 * @param lp unused
 * @return Zero-sum vector or matrix which is one larger in each dimension
 */
template <typename T, typename Lp, require_t<is_rev_matrix<T>>* = nullptr>
inline auto sum_to_zero_constrain(T&& y, Lp& lp) {
  return sum_to_zero_constrain(std::forward<T>(y));
}

}  // namespace math
}  // namespace stan
#endif
