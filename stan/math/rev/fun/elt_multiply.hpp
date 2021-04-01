#ifndef STAN_MATH_REV_FUN_ELT_MULTIPLY_HPP
#define STAN_MATH_REV_FUN_ELT_MULTIPLY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/multiply.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise multiplication of the specified
 * matrices.
 *
 * @tparam Mat1 type of the first matrix or expression
 * @tparam Mat2 type of the second matrix or expression
 *
 * @param m1 First matrix or expression
 * @param m2 Second matrix or expression
 * @return Elementwise product of matrices.
 */
template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2>* = nullptr,
          require_any_rev_matrix_t<Mat1, Mat2>* = nullptr>
auto elt_multiply(const Mat1& m1, const Mat2& m2) {
  check_matching_dims("elt_multiply", "m1", m1, "m2", m2);
  using inner_ret_type = decltype(value_of(m1).cwiseProduct(value_of(m2)));
  using ret_type = return_var_matrix_t<inner_ret_type, Mat1, Mat2>;
  if (!is_constant<Mat1>::value && !is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_m1 = m1;
    arena_t<promote_scalar_t<var, Mat2>> arena_m2 = m2;
    arena_t<ret_type> ret(arena_m1.val().cwiseProduct(arena_m2.val()));
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      for (Eigen::Index j = 0; j < arena_m2.cols(); ++j) {
        for (Eigen::Index i = 0; i < arena_m2.rows(); ++i) {
          const auto ret_adj = ret.adj().coeffRef(i, j);
          arena_m1.adj().coeffRef(i, j) += arena_m2.val().coeff(i, j) * ret_adj;
          arena_m2.adj().coeffRef(i, j) += arena_m1.val().coeff(i, j) * ret_adj;
        }
      }
    });
    return ret_type(ret);
  } else if (!is_constant<Mat1>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_m1 = m1;
    arena_t<promote_scalar_t<double, Mat2>> arena_m2 = value_of(m2);
    arena_t<ret_type> ret(arena_m1.val().cwiseProduct(arena_m2));
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      arena_m1.adj().array() += arena_m2.array() * ret.adj().array();
    });
    return ret_type(ret);
  } else if (!is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<double, Mat1>> arena_m1 = value_of(m1);
    arena_t<promote_scalar_t<var, Mat2>> arena_m2 = m2;
    arena_t<ret_type> ret(arena_m1.cwiseProduct(arena_m2.val()));
    reverse_pass_callback([ret, arena_m2, arena_m1]() mutable {
      arena_m2.adj().array() += arena_m1.array() * ret.adj().array();
    });
    return ret_type(ret);
  }
}

}  // namespace math
}  // namespace stan

#endif
