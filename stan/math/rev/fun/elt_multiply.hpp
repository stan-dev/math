#ifndef STAN_MATH_REV_FUN_ELT_MULTIPLY_HPP
#define STAN_MATH_REV_FUN_ELT_MULTIPLY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/eval.hpp>
#include <stan/math/rev/core.hpp>

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
  const auto& m1_ref = to_ref(m1);
  const auto& m2_ref = to_ref(m2);
  check_matching_dims("elt_multiply", "m1", m1_ref, "m2", m2_ref);
  using inner_ret_type
      = decltype(value_of(m1_ref).cwiseProduct(value_of(m2_ref)));
  using ret_type = promote_var_matrix_t<inner_ret_type, Mat1, Mat2>;
  if (!is_constant<Mat1>::value && !is_constant<Mat2>::value) {
    arena_t<Mat1> arena_m1 = m1_ref;
    arena_t<Mat2> arena_m2 = m2_ref;
    arena_t<ret_type> ret(value_of(arena_m1).cwiseProduct(value_of(arena_m2)));
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      using var_m1 = arena_t<promote_scalar_t<var, Mat1>>;
      using var_m2 = arena_t<promote_scalar_t<var, Mat2>>;
      if (is_var_matrix<Mat1>::value || is_var_matrix<Mat2>::value) {
        forward_as<var_m1>(arena_m1).adj()
            += value_of(arena_m2).cwiseProduct(ret.adj());
        forward_as<var_m2>(arena_m2).adj()
            += value_of(arena_m1).cwiseProduct(ret.adj());
      } else {
        for (Eigen::Index i = 0; i < arena_m2.size(); ++i) {
          forward_as<var_m1>(arena_m1).coeffRef(i).adj()
              += forward_as<var_m2>(arena_m2).coeffRef(i).val()
                 * ret.coeffRef(i).adj();
          forward_as<var_m2>(arena_m2).coeffRef(i).adj()
              += forward_as<var_m1>(arena_m1).coeffRef(i).val()
                 * ret.coeffRef(i).adj();
        }
      }
    });
    return ret_type(ret);
  } else if (!is_constant<Mat1>::value) {
    arena_t<Mat1> arena_m1 = m1_ref;
    auto arena_m2 = to_arena(value_of(m2_ref));
    arena_t<ret_type> ret(value_of(arena_m1).cwiseProduct(arena_m2));
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      using var_m1 = arena_t<promote_scalar_t<var, Mat1>>;
      forward_as<var_m1>(arena_m1).adj().array()
          += arena_m2.array() * ret.adj().array();
    });
    return ret_type(ret);
  } else if (!is_constant<Mat2>::value) {
    arena_t<Mat2> arena_m2 = m2_ref;
    auto arena_m1 = to_arena(value_of(m1_ref));
    arena_t<ret_type> ret(arena_m1.cwiseProduct(value_of(arena_m2)));
    reverse_pass_callback([ret, arena_m2, arena_m1]() mutable {
      using var_m2 = arena_t<promote_scalar_t<var, Mat2>>;
      forward_as<var_m2>(arena_m2).adj().array()
          += arena_m1.array() * ret.adj().array();
    });
    return ret_type(ret);
  }
}
}  // namespace math
}  // namespace stan

#endif
