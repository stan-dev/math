#ifndef STAN_MATH_REV_FUN_FROM_VAR_VALUE_HPP
#define STAN_MATH_REV_FUN_FROM_VAR_VALUE_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core/arena_matrix.hpp>
#include <stan/math/rev/core/reverse_pass_callback.hpp>

namespace stan {
namespace math {

/**
 * Converts `var_value` into an Eigen Matrix. Adjoint is propagated back to
 * argument in the reverse pass.
 *
 * @tparam T type of the input
 * @param a matrix to convert
 */
template <typename T, require_var_matrix_t<T>* = nullptr>
Eigen::Matrix<var, T::RowsAtCompileTime, T::ColsAtCompileTime> from_var_value(
    const T& a) {
  arena_matrix<Eigen::Matrix<var, T::RowsAtCompileTime, T::ColsAtCompileTime>>
      res(a.val());
  reverse_pass_callback([res, a]() mutable { a.vi_->adj_ += res.adj(); });
  return res;
}

/**
 * This is a no-op for Eigen containers of vars, scalars or prim types.
 *
 * @tparam T type of the input
 * @param a matrix to convert
 */
template <
    typename T,
    require_any_t<
        conjunction<is_eigen<T>, is_var<scalar_type_t<T>>>,
        std::is_same<std::decay_t<T>, var>,
        bool_constant<!std::is_same<scalar_type_t<T>, var>::value>>* = nullptr>
T from_var_value(T&& a) {
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
auto from_var_value(const std::vector<T>& a) {
  std::vector<decltype(from_var_value(std::declval<T>()))> out;
  out.reserve(a.size());
  for (size_t i = 0; i < a.size(); ++i) {
    out.emplace_back(from_var_value(a[i]));
  }
  return out;
}

}  // namespace math
}  // namespace stan

#endif
