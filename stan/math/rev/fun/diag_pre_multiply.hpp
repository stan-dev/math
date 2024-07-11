#ifndef STAN_MATH_REV_FUN_DIAG_PRE_MULTIPLY_HPP
#define STAN_MATH_REV_FUN_DIAG_PRE_MULTIPLY_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/rev/core.hpp>

namespace stan {
namespace math {

/**
 * Return the product of the diagonal matrix formed from the vector
 * or row_vector and a matrix.
 *
 * @tparam T1 type of the vector/row_vector
 * @tparam T2 type of the matrix
 * @param m1 input vector/row_vector
 * @param m2 input matrix
 *
 * @return product of the diagonal matrix formed from the
 * vector or row_vector and a matrix.
 */
template <typename T1, typename T2, require_vector_t<T1>* = nullptr,
          require_matrix_t<T2>* = nullptr,
          require_any_st_var<T1, T2>* = nullptr>
inline auto diag_pre_multiply(T1&& m1, T2&& m2) {
  check_size_match("diag_pre_multiply", "m1.size()", m1.size(), "m2.rows()",
                   m2.rows());
  using inner_ret_type = decltype(value_of(m1).asDiagonal() * value_of(m2));
  using ret_type = return_var_matrix_t<inner_ret_type, T1, T2>;
  arena_t<T1> arena_m1 = std::forward<T1>(m1);
  arena_t<T2> arena_m2 = std::forward<T2>(m2);
  if constexpr (!is_constant_v<T1> && !is_constant_v<T2>) {
    arena_t<ret_type> ret(arena_m1.val().asDiagonal() * arena_m2.val());
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      arena_m1.adj() += arena_m2.val().cwiseProduct(ret.adj()).rowwise().sum();
      arena_m2.adj() += arena_m1.val().asDiagonal() * ret.adj();
    });
    return ret;
  } else if constexpr (!is_constant_v<T1>) {
    arena_t<ret_type> ret(arena_m1.val().asDiagonal() * arena_m2);
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      arena_m1.adj() += arena_m2.val().cwiseProduct(ret.adj()).rowwise().sum();
    });
    return ret;
  } else if constexpr (!is_constant_v<T2>) {
    arena_t<ret_type> ret(arena_m1.asDiagonal() * arena_m2.val());
    reverse_pass_callback([ret, arena_m1, arena_m2]() mutable {
      arena_m2.adj() += arena_m1.val().asDiagonal() * ret.adj();
    });
    return ret;
  }
}

}  // namespace math
}  // namespace stan

#endif
