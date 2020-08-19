#ifndef STAN_MATH_REV_FUN_COLUMNS_DOT_PRODUCT_HPP
#define STAN_MATH_REV_FUN_COLUMNS_DOT_PRODUCT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/functor/arena_matrix.hpp>
#include <stan/math/prim/meta.hpp>

#include <type_traits>

namespace stan {
namespace math {

template <typename Mat1, typename Mat2,
          require_all_eigen_t<Mat1, Mat2>* = nullptr,
          require_any_eigen_vt<is_var, Mat1, Mat2>* = nullptr>
inline Eigen::Matrix<return_type_t<Mat1, Mat2>, 1, Mat1::ColsAtCompileTime>
columns_dot_product(const Mat1& v1, const Mat2& v2) {
  check_matching_sizes("columns_dot_product", "v1", v1, "v2", v2);

  ref_type_t<Mat1> v1_ref = v1;
  ref_type_t<Mat2> v2_ref = v2;

  arena_matrix<promote_scalar_t<double, Mat1>> arena_v1_val = value_of(v1_ref);
  arena_matrix<promote_scalar_t<double, Mat2>> arena_v2_val = value_of(v2_ref);

  arena_matrix<promote_scalar_t<var, Mat1>> arena_v1;
  arena_matrix<promote_scalar_t<var, Mat2>> arena_v2;

  if (!is_constant<Mat1>::value) {
    arena_v1 = v1_ref;
  }

  if (!is_constant<Mat2>::value) {
    arena_v2 = v2_ref;
  }

  Eigen::RowVectorXd out_val(arena_v1_val.cols());
  for (size_t m = 0; m < arena_v1_val.cols(); ++m)
    out_val.coeffRef(m) = arena_v1_val.col(m).dot(arena_v2_val.col(m));

  arena_matrix<Eigen::Matrix<var, 1, Eigen::Dynamic>> out = out_val;

  reverse_pass_callback(
      [out, arena_v1, arena_v2, arena_v1_val, arena_v2_val]() mutable {
        Eigen::RowVectorXd adj = out.adj();

        if (!is_constant<Mat1>::value)
          arena_v1.adj() += arena_v2_val * adj.asDiagonal();

        if (!is_constant<Mat2>::value)
          arena_v2.adj() += arena_v1_val * adj.asDiagonal();
      });

  return out;
}

}  // namespace math
}  // namespace stan
#endif
