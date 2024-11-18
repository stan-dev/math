#ifndef STAN_MATH_REV_FUN_COLUMNS_DOT_PRODUCT_HPP
#define STAN_MATH_REV_FUN_COLUMNS_DOT_PRODUCT_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/rev/fun/dot_product.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>

#include <type_traits>

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
          require_any_eigen_vt<is_var, Mat1, Mat2>* = nullptr>
inline Eigen::Matrix<return_type_t<Mat1, Mat2>, 1, Mat1::ColsAtCompileTime>
columns_dot_product(const Mat1& v1, const Mat2& v2) {
  check_matching_sizes("dot_product", "v1", v1, "v2", v2);
  Eigen::Matrix<var, 1, Mat1::ColsAtCompileTime> ret(1, v1.cols());
  for (size_type j = 0; j < v1.cols(); ++j) {
    ret.coeffRef(j) = dot_product(v1.col(j), v2.col(j));
  }
  return ret;
}

/**
 * Returns the dot product of columns of the specified matrices.
 *
 * This overload is used when at least one of Mat1 and Mat2 is
 * a `var_value<T>` where `T` inherits from `EigenBase`. The other
 * argument can be another `var_value` or a type that inherits from
 * `EigenBase`.
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
          require_all_matrix_t<Mat1, Mat2>* = nullptr,
          require_any_var_matrix_t<Mat1, Mat2>* = nullptr>
inline auto columns_dot_product(const Mat1& v1, const Mat2& v2) {
  check_matching_sizes("columns_dot_product", "v1", v1, "v2", v2);
  using inner_return_t = decltype(
      (value_of(v1).array() * value_of(v2).array()).colwise().sum().matrix());
  using return_t = return_var_matrix_t<inner_return_t, Mat1, Mat2>;

  if (!is_constant<Mat1>::value && !is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_v1 = v1;
    arena_t<promote_scalar_t<var, Mat2>> arena_v2 = v2;

    return_t res
        = (arena_v1.val().array() * arena_v2.val().array()).colwise().sum();

    reverse_pass_callback([arena_v1, arena_v2, res]() mutable {
      if (is_var_matrix<Mat1>::value) {
        arena_v1.adj().noalias() += arena_v2.val() * res.adj().asDiagonal();
      } else {
        arena_v1.adj() += arena_v2.val() * res.adj().asDiagonal();
      }
      if (is_var_matrix<Mat2>::value) {
        arena_v2.adj().noalias() += arena_v1.val() * res.adj().asDiagonal();
      } else {
        arena_v2.adj() += arena_v1.val() * res.adj().asDiagonal();
      }
    });

    return res;
  } else if (!is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<double, Mat1>> arena_v1 = value_of(v1);
    arena_t<promote_scalar_t<var, Mat2>> arena_v2 = v2;

    return_t res = (arena_v1.array() * arena_v2.val().array()).colwise().sum();

    reverse_pass_callback([arena_v1, arena_v2, res]() mutable {
      if (is_var_matrix<Mat2>::value) {
        arena_v2.adj().noalias() += arena_v1 * res.adj().asDiagonal();
      } else {
        arena_v2.adj() += arena_v1 * res.adj().asDiagonal();
      }
    });

    return res;
  } else {
    arena_t<promote_scalar_t<var, Mat1>> arena_v1 = v1;
    arena_t<promote_scalar_t<double, Mat2>> arena_v2 = value_of(v2);

    return_t res = (arena_v1.val().array() * arena_v2.array()).colwise().sum();

    reverse_pass_callback([arena_v1, arena_v2, res]() mutable {
      if (is_var_matrix<Mat2>::value) {
        arena_v1.adj().noalias() += arena_v2 * res.adj().asDiagonal();
      } else {
        arena_v1.adj() += arena_v2 * res.adj().asDiagonal();
      }
    });

    return res;
  }
}

}  // namespace math
}  // namespace stan
#endif
