#ifndef STAN_MATH_REV_FUN_COLUMNS_DOT_PRODUCT_HPP
#define STAN_MATH_REV_FUN_COLUMNS_DOT_PRODUCT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/prim/meta.hpp>

#include <type_traits>

namespace stan {
namespace math {

/**
 * Returns the dot product of columns of the specified matrices.
 *
 * @tparam Mat1 type of the first matrix
 * @tparam Mat2 type of the second matrix
 *
 * @param v1 Matrix of first vectors.
 * @param v2 Matrix of second vectors.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2>* = nullptr>
inline auto columns_dot_product(Mat1&& v1, Mat2&& v2) {
  check_matching_sizes("columns_dot_product", "v1", v1, "v2", v2);
  using inner_return_t = decltype(
      (value_of(v1).array() * value_of(v2).array()).colwise().sum().matrix());
  using return_t = return_var_matrix_t<inner_return_t, Mat1, Mat2>;
  arena_t<Mat1> arena_v1 = std::forward<Mat1>(v1);
  arena_t<Mat2> arena_v2 = std::forward<Mat2>(v2);
  arena_t<return_t> res
      = (value_of(arena_v1).array() * value_of(arena_v2).array())
            .colwise()
            .sum();
  reverse_pass_callback([arena_v1, arena_v2, res]() mutable {
    if constexpr (is_autodiffable_v<Mat1>) {
      if constexpr (is_var_matrix<Mat1>::value) {
        arena_v1.adj().noalias() += value_of(arena_v2) * res.adj().asDiagonal();
      } else {
        arena_v1.adj() += value_of(arena_v2) * res.adj().asDiagonal();
      }
    }
    if constexpr (is_autodiffable_v<Mat2>) {
      if constexpr (is_var_matrix<Mat2>::value) {
        arena_v2.adj().noalias() += value_of(arena_v1) * res.adj().asDiagonal();
      } else {
        arena_v2.adj() += value_of(arena_v1) * res.adj().asDiagonal();
      }
    }
  });
  return res;
}

}  // namespace math
}  // namespace stan
#endif
