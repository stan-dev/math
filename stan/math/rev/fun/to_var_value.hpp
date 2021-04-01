#ifndef STAN_MATH_REV_FUN_TO_VAR_VALUE_HPP
#define STAN_MATH_REV_FUN_TO_VAR_VALUE_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>

namespace stan {
namespace math {

/**
 * Converts an Eigen matrix (or vector or row_vector) or expression of `var`s
 * into `var_value`. Adjoint is propagated back to argument in the reverse
 * pass.
 *
 * @tparam T type of the input
 * @param a matrix to convert
 */
template <typename T, require_eigen_vt<is_var, T>* = nullptr>
var_value<Eigen::Matrix<double, T::RowsAtCompileTime, T::ColsAtCompileTime>>
to_var_value(const T& a) {
  arena_matrix<plain_type_t<T>> a_arena = a;
  var_value<promote_scalar_t<double, T>> res(a_arena.val());
  reverse_pass_callback(
      [res, a_arena]() mutable { a_arena.adj() += res.adj(); });
  return res;
}

/**
 * This is a no-op for var_values.
 *
 * @tparam T type of the input
 * @param a matrix to convert
 */
template <typename T, require_var_t<T>* = nullptr>
T to_var_value(T&& a) {
  return std::forward<T>(a);
}

/**
 * Convert the elements of the `std::vector` input to `var_value` types
 * if possible
 *
 * @tparam T type of elemnts of the input vector
 * @param a std::vector of elements to convert
 */
template <typename T>
auto to_var_value(const std::vector<T>& a) {
  std::vector<decltype(to_var_value(std::declval<T>()))> out;
  out.reserve(a.size());
  for (size_t i = 0; i < a.size(); ++i) {
    out.emplace_back(to_var_value(a[i]));
  }
  return out;
}

}  // namespace math
}  // namespace stan

#endif
