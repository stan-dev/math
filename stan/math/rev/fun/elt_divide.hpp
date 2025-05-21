#ifndef STAN_MATH_REV_FUN_ELT_DIVIDE_HPP
#define STAN_MATH_REV_FUN_ELT_DIVIDE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/prim/fun/elt_divide.hpp>

namespace stan {
namespace math {

/**
 * Return the elementwise division of the specified
 * matrices.
 *
 * @tparam Mat1 type of the first matrix
 * @tparam Mat2 type of the second matrix
 *
 * @param m1 First matrix
 * @param m2 Second matrix
 * @return Elementwise division of matrices.
 */
template <typename Mat1, typename Mat2,
          require_all_matrix_t<Mat1, Mat2>* = nullptr,
          require_any_rev_matrix_t<Mat1, Mat2>* = nullptr>
auto elt_divide(const Mat1& m1, const Mat2& m2) {
  check_matching_dims("elt_divide", "m1", m1, "m2", m2);
  using inner_ret_type
      = decltype((value_of(m1).array() / value_of(m2).array()).matrix());
  using ret_type = return_var_matrix_t<inner_ret_type, Mat1, Mat2>;
  if (!is_constant<Mat1>::value && !is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_m1 = m1;
    arena_t<promote_scalar_t<var, Mat2>> arena_m2 = m2;
    arena_t<ret_type> ret(arena_m1.val().array() / arena_m2.val().array());
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      for (Eigen::Index j = 0; j < arena_m2.cols(); ++j) {
        for (Eigen::Index i = 0; i < arena_m2.rows(); ++i) {
          const auto ret_div
              = ret.adj().coeff(i, j) / arena_m2.val().coeff(i, j);
          arena_m1.adj().coeffRef(i, j) += ret_div;
          arena_m2.adj().coeffRef(i, j) -= ret.val().coeff(i, j) * ret_div;
        }
      }
    });
    return ret_type(ret);
  } else if (!is_constant<Mat1>::value) {
    arena_t<promote_scalar_t<var, Mat1>> arena_m1 = m1;
    arena_t<promote_scalar_t<double, Mat2>> arena_m2 = value_of(m2);
    arena_t<ret_type> ret(arena_m1.val().array() / arena_m2.array());
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      arena_m1.adj().array() += ret.adj().array() / arena_m2.array();
    });
    return ret_type(ret);
  } else if (!is_constant<Mat2>::value) {
    arena_t<promote_scalar_t<double, Mat1>> arena_m1 = value_of(m1);
    arena_t<promote_scalar_t<var, Mat2>> arena_m2 = m2;
    arena_t<ret_type> ret(arena_m1.array() / arena_m2.val().array());
    reverse_pass_callback([ret, arena_m2, arena_m1]() mutable {
      arena_m2.adj().array()
          -= ret.val().array() * ret.adj().array() / arena_m2.val().array();
    });
    return ret_type(ret);
  }
}

/**
 * Return the elementwise division of the specified scalar
 * by the specified matrix.
 *
 * @tparam Scal type of the scalar
 * @tparam Mat type of the matrix
 *
 * @param s scalar
 * @param m matrix or expression
 * @return Elementwise division of a scalar by matrix.
 */
template <typename Scal, typename Mat, require_stan_scalar_t<Scal>* = nullptr,
          require_var_matrix_t<Mat>* = nullptr>
auto elt_divide(Scal s, const Mat& m) {
  plain_type_t<Mat> res = value_of(s) / m.val().array();

  reverse_pass_callback([m, s, res]() mutable {
    m.adj().array() -= res.val().array() * res.adj().array() / m.val().array();
    if (!is_constant<Scal>::value)
      forward_as<var>(s).adj() += (res.adj().array() / m.val().array()).sum();
  });

  return res;
}

}  // namespace math
}  // namespace stan

#endif
