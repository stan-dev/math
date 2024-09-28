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
    for (int i = 0; i < N; ++i) {
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
template <typename T, require_rev_col_vector_t<T>* = nullptr>
inline auto sum_to_zero_constrain(T&& y, scalar_type_t<T>& lp) {
  return sum_to_zero_constrain(std::forward<T>(y));
}

}  // namespace math
}  // namespace stan
#endif
