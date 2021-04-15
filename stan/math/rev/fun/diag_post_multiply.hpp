#ifndef STAN_MATH_REV_FUN_DIAG_POST_MULTIPLY_HPP
#define STAN_MATH_REV_FUN_DIAG_POST_MULTIPLY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Return the product of the matrix and a diagonal matrix formed from the vector
 * or row_vector.
 *
 * @tparam T1 type of the matrix
 * @tparam T2 type of the vector/row_vector
 * @param m1 input matrix
 * @param m2 input vector/row_vector
 *
 * @return product of the matrix and the diagonal matrix formed from the
 * vector or row_vector.
 */
template <typename T1, typename T2, require_matrix_t<T1>* = nullptr,
          require_vector_t<T2>* = nullptr,
          require_any_st_var<T1, T2>* = nullptr>
auto diag_post_multiply(const T1& m1, const T2& m2) {
  check_size_match("diag_post_multiply", "m2.size()", m2.size(), "m1.cols()",
                   m1.cols());
  using inner_ret_type = decltype(value_of(m1) * value_of(m2).asDiagonal());
  using ret_type = return_var_matrix_t<inner_ret_type, T1, T2>;

  if (!is_constant<T1>::value && !is_constant<T2>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_m1 = m1;
    arena_t<promote_scalar_t<var, T2>> arena_m2 = m2;
    arena_t<ret_type> ret(arena_m1.val() * arena_m2.val().asDiagonal());
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      arena_m2.adj() += arena_m1.val().cwiseProduct(ret.adj()).colwise().sum();
      arena_m1.adj() += ret.adj() * arena_m2.val().asDiagonal();
    });
    return ret_type(ret);
  } else if (!is_constant<T1>::value) {
    arena_t<promote_scalar_t<var, T1>> arena_m1 = m1;
    arena_t<promote_scalar_t<double, T2>> arena_m2 = value_of(m2);
    arena_t<ret_type> ret(arena_m1.val() * arena_m2.asDiagonal());
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      arena_m1.adj() += ret.adj() * arena_m2.val().asDiagonal();
    });
    return ret_type(ret);
  } else if (!is_constant<T2>::value) {
    arena_t<promote_scalar_t<double, T1>> arena_m1 = value_of(m1);
    arena_t<promote_scalar_t<var, T2>> arena_m2 = m2;
    arena_t<ret_type> ret(arena_m1 * arena_m2.val().asDiagonal());
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      arena_m2.adj() += arena_m1.val().cwiseProduct(ret.adj()).colwise().sum();
    });
    return ret_type(ret);
  }
}

}  // namespace math
}  // namespace stan

#endif
