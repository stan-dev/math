#ifndef STAN_MATH_REV_FUN_COLUMNS_DOT_PRODUCT_HPP
#define STAN_MATH_REV_FUN_COLUMNS_DOT_PRODUCT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/typedefs.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
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

  arena_matrix<promote_scalar_t<var, Mat1>> arena_v1 = to_arena_if<!is_constant<Mat1>::value>(v1_ref);
  arena_matrix<promote_scalar_t<var, Mat2>> arena_v2 = to_arena_if<!is_constant<Mat2>::value>(v2_ref);

  Eigen::RowVectorXd out_val = (arena_v1_val.cwiseProduct(arena_v2_val)).colwise().sum();
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

/**
 * Deduce whether the return type should be a Matrix<var> or var<Matrix>
 * There rule here is that for the supplied return type we want, if
 * any of the types are Eigen matrices then the return type is just ReturnType.
 * But if none of them are then the return
 */
template <typename ReturnType, typename... Args>
using return_var_matrix_t = std::conditional_t<disjunction<is_eigen<std::decay_t<Args>>...>::value,
  ReturnType, var_matrix_converter_t<ReturnType>>;


template <typename T>
using is_var_eigen = bool_constant<//(is_eigen<T>::value && is_var<value_type_t<T>>::value) ||
 (is_var<T>::value && is_eigen<value_type_t<T>>::value)>;

template <typename Mat1, typename Mat2,
          require_any_t<is_var_eigen<Mat1>, is_var_eigen<Mat2>>* = nullptr>
inline var_value<Eigen::Matrix<double, 1, Mat1::RowsAtCompileTime>> columns_dot_product(const Mat1& v1, const Mat2& v2) {
  check_matching_sizes("columns_dot_product", "v1", v1, "v2", v2);

  arena_matrix<var_matrix_converter_t<Mat1>> arena_v1 = to_arena(v1);
  arena_matrix<var_matrix_converter_t<Mat2>> arena_v2 = to_arena(v2);

  Eigen::Matrix<double, 1, Mat1::ColsAtCompileTime> out_val = (arena_v1.val().cwiseProduct(arena_v2.val())).colwise().sum();
  arena_matrix<var_value<Eigen::Matrix<double, 1, Mat1::RowsAtCompileTime>>> out(out_val);

  reverse_pass_callback(
      [out, arena_v1, arena_v2]() mutable {
        Eigen::RowVectorXd adj = out.adj();

        if (!is_constant<Mat1>::value)
          arena_v1.adj() += arena_v2.val() * adj.asDiagonal();

        if (!is_constant<Mat2>::value)
          arena_v2.adj() += arena_v1.val() * adj.asDiagonal();
      });

  return out;
}

}  // namespace math
}  // namespace stan
#endif
