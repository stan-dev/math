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
inline auto rows_dot_product(const Mat1& v1, const Mat2& v2) {
  using ret_type = promote_scalar_t<var, Eigen::VectorXd>;

  check_matching_sizes("rows_dot_product", "v1", v1, "v2", v2);

  if (!is_constant<Mat1>::value && !is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_v1 = v1;
    arena_t<promote_scalar_t<var, Mat2>> arena_v2 = v2;

    auto arena_v1_val = to_arena(value_of(arena_v1));
    auto arena_v2_val = to_arena(value_of(arena_v2));

    arena_t<ret_type> res(arena_v1_val.rows());
    for (size_t m = 0; m < arena_v1_val.rows(); ++m)
      res.coeffRef(m) = arena_v1_val.row(m).dot(arena_v2_val.row(m));

    reverse_pass_callback(
        [res, arena_v1, arena_v2, arena_v1_val, arena_v2_val]() mutable {
          for (size_t j = 0; j < arena_v1_val.cols(); ++j) {
            for (size_t i = 0; i < arena_v1_val.rows(); ++i) {
              arena_v1.coeffRef(i, j).adj()
                  += arena_v2_val.coeff(i, j) * res.coeff(i).adj();
              arena_v2.coeffRef(i, j).adj()
                  += arena_v1_val.coeff(i, j) * res.coeff(i).adj();
            }
          }
        });

    return ret_type(res);
  } else if (!is_constant<Mat1>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_v1 = v1;

    auto arena_v2_val = to_arena(value_of(v2));

    arena_t<ret_type> res(arena_v2_val.rows());
    for (size_t m = 0; m < arena_v2_val.rows(); ++m)
      res.coeffRef(m) = arena_v1.row(m).val().dot(arena_v2_val.row(m));

    reverse_pass_callback([res, arena_v1, arena_v2_val]() mutable {
      for (size_t j = 0; j < arena_v2_val.cols(); ++j) {
        for (size_t i = 0; i < arena_v2_val.rows(); ++i) {
          arena_v1.coeffRef(i, j).adj()
              += arena_v2_val.coeff(i, j) * res.coeff(i).adj();
        }
      }
    });

    return ret_type(res);
  } else {
    arena_t<promote_scalar_t<var, Mat2>> arena_v2 = v2;

    auto arena_v1_val = to_arena(value_of(v1));

    arena_t<ret_type> res(arena_v1_val.rows());
    for (size_t m = 0; m < arena_v1_val.rows(); ++m)
      res.coeffRef(m) = arena_v1_val.row(m).dot(arena_v2.row(m).val());

    reverse_pass_callback([res, arena_v2, arena_v1_val]() mutable {
      for (size_t j = 0; j < arena_v1_val.cols(); ++j) {
        for (size_t i = 0; i < arena_v1_val.rows(); ++i) {
          arena_v2.coeffRef(i, j).adj()
              += arena_v1_val.coeff(i, j) * res.coeff(i).adj();
        }
      }
    });

    return ret_type(res);
  }
}

}  // namespace math
}  // namespace stan
#endif
