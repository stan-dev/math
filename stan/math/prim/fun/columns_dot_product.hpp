#ifndef STAN_MATH_PRIM_FUN_COLUMNS_DOT_PRODUCT_HPP
#define STAN_MATH_PRIM_FUN_COLUMNS_DOT_PRODUCT_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns the dot product of columns of the specified matrices.
 *
 * @tparam Mat1 type of the first matrix (must be derived from \c
 * Eigen::MatrixBase)
 * @tparam Mat2 type of the second matrix (must be derived from \c
 * Eigen::MatrixBase)
 *
 * @param v1 Matrix of first vectors.
 * @param v2 Matrix of second vectors.
 * @return Dot product of the vectors.
 * @throw std::domain_error If the vectors are not the same
 * size or if they are both not vector dimensioned.
 */
template <typename Mat1, typename Mat2,
          require_all_eigen_t<Mat1, Mat2>* = nullptr,
          require_all_not_eigen_vt<is_var, Mat1, Mat2>* = nullptr>
inline auto columns_dot_product(Mat1&& v1, Mat2&& v2) {
  check_matching_dims("columns_dot_product", "v1", v1, "v2", v2);
  return make_holder(
      [](auto&& v1_, auto&& v2_) {
        return std::forward<decltype(v1_)>(v1_)
            .cwiseProduct(std::forward<decltype(v2_)>(v2_))
            .colwise()
            .sum();
      },
      std::forward<Mat1>(v1), std::forward<Mat2>(v2));
}

}  // namespace math
}  // namespace stan

#endif
