#ifndef STAN_MATH_REV_FUN_COLUMNS_DOT_PRODUCT_HPP
#define STAN_MATH_REV_FUN_COLUMNS_DOT_PRODUCT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/functor/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/prim/meta.hpp>

#include <type_traits>

namespace stan {
namespace math {

template <typename Mat1, typename Mat2,
          require_all_eigen_t<Mat1, Mat2>* = nullptr,
          require_any_eigen_vt<is_var, Mat1, Mat2>* = nullptr>
inline auto columns_dot_product(const Mat1& v1, const Mat2& v2) {
  using ret_type = promote_scalar_t<var, Eigen::RowVectorXd>;
  using mat1_var_type = promote_scalar_t<var, Mat1>;
  using mat2_var_type = promote_scalar_t<var, Mat2>;

  check_matching_sizes("columns_dot_product", "v1", v1, "v2", v2);

  auto arena_v1 = to_arena(v1);
  auto arena_v2 = to_arena(v2);

  arena_t<ret_type> res
      = (arena_v1.val().array() * arena_v2.val().array()).colwise().sum();

  reverse_pass_callback([res, arena_v1, arena_v2]() mutable {
    for (size_t j = 0; j < arena_v1.cols(); ++j) {
      for (size_t i = 0; i < arena_v1.rows(); ++i) {
        if (!is_constant<Mat1>::value)
          forward_as<var>(arena_v1.coeffRef(i, j)).adj()
              += value_of(arena_v2.coeff(i, j)) * res.coeff(j).adj();
        if (!is_constant<Mat2>::value)
          forward_as<var>(arena_v2.coeffRef(i, j)).adj()
              += value_of(arena_v1.coeff(i, j)) * res.coeff(j).adj();
      }
    }
  });

  return ret_type(res);
}

}  // namespace math
}  // namespace stan
#endif
