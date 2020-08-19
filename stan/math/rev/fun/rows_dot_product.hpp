#ifndef STAN_MATH_REV_FUN_ROWS_DOT_PRODUCT_HPP
#define STAN_MATH_REV_FUN_ROWS_DOT_PRODUCT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <type_traits>

namespace stan {
namespace math {

template <typename Mat1, typename Mat2,
          require_all_eigen_t<Mat1, Mat2>* = nullptr,
          require_any_eigen_vt<is_var, Mat1, Mat2>* = nullptr>
inline Eigen::Matrix<var, Mat1::RowsAtCompileTime, 1> rows_dot_product(
    const Mat1& v1, const Mat2& v2) {
  check_matching_sizes("rows_dot_product", "v1", v1, "v2", v2);

  ref_type_t<Mat1> v1_ref = v1;
  ref_type_t<Mat2> v2_ref = v2;
  
  using Mat1_double = plain_type_t<decltype(value_of(v1_ref))>;
  using Mat1_var = promote_scalar_t<var, Mat1_double>;

  using Mat2_double = plain_type_t<decltype(value_of(v2_ref))>;
  using Mat2_var = promote_scalar_t<var, Mat2_double>;

  arena_matrix<Mat1_double> arena_v1_val = value_of(v1_ref);
  arena_matrix<Mat2_double> arena_v2_val = value_of(v2_ref);

  arena_matrix<Mat1_var> arena_v1;
  arena_matrix<Mat2_var> arena_v2;

  if(!is_constant<Mat1>::value) {
    arena_v1 = v1_ref;
  }

  if(!is_constant<Mat2>::value) {
    arena_v2 = v2_ref;
  }

  Eigen::VectorXd out_val(arena_v1_val.rows());
  for(size_t m = 0; m < arena_v1_val.rows(); ++m)
    out_val.coeffRef(m) = arena_v1_val.row(m).dot(arena_v2_val.row(m));

  arena_matrix<Eigen::Matrix<var, Eigen::Dynamic, 1>> out = out_val;

  reverse_pass_callback([out, arena_v1, arena_v2,
			 arena_v1_val, arena_v2_val]() mutable {
    Eigen::VectorXd adj = out.adj();

    if(!is_constant<Mat1>::value)
      arena_v1.adj() += adj.asDiagonal() * arena_v2_val;

    if(!is_constant<Mat2>::value)
      arena_v2.adj() += adj.asDiagonal() * arena_v1_val;
  });

  return out;
}

}  // namespace math
}  // namespace stan
#endif
