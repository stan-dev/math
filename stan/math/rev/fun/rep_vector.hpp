#ifndef STAN_MATH_REV_FUN_REP_VECTOR_HPP
#define STAN_MATH_REV_FUN_REP_VECTOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/err.hpp>

namespace stan {
namespace math {

/**
 * Overload for `var_value<Vector>`.
 * @tparam T_ret The user supplied return type.
 * @tparam T A double or var type
 * @param x The type to be propogated through the new vector.
 * @param n The size of the new vector.
 */
template <typename T_ret, typename T, require_var_matrix_t<T_ret>* = nullptr,
          require_eigen_col_vector_t<value_type_t<T_ret>>* = nullptr,
          require_stan_scalar_t<T>* = nullptr>
inline auto rep_vector(const T& x, int n) {
  check_nonnegative("rep_vector", "n", n);
  return make_callback_var(value_type_t<T_ret>::Constant(n, value_of(x)),
                           [x](auto& vi) mutable {
                             if (is_var<T>::value) {
                               forward_as<var>(x).adj() += vi.adj().sum();
                             }
                           });
}

}  // namespace math
}  // namespace stan

#endif
