#ifndef STAN_MATH_REV_MAT_FUN_TO_VAR_HPP
#define STAN_MATH_REV_MAT_FUN_TO_VAR_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/typedefs.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/fun/typedefs.hpp>
#include <stan/math/rev/scal/fun/to_var.hpp>

namespace stan {
namespace math {

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] x A Matrix with scalars
 * @return A Matrix with automatic differentiation variables
 */
template <typename T, enable_if_eigen<T>...,
          enable_if_var<scalar_type_decay_t<T>>...>
inline auto&& to_var(T&& x) {
  return std::forward<T>(x);
}

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] x A Matrix with scalars
 * @return A Matrix with automatic differentiation variables
 */
inline matrix_v to_var(const matrix_d& x) {
  matrix_v m_v = x;
  return m_v;
}

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] x A Vector of scalars
 * @return A Vector of automatic differentiation variables with
 *   values of x
 */
inline vector_v to_var(const vector_d& x) {
  vector_v v_v = x;
  return v_v;
}

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] x A row vector of scalars
 * @return A row vector of automatic differentation variables with
 *   values of x.
 */
inline row_vector_v to_var(const row_vector_d& x) {
  row_vector_v rv_v = x;
  return rv_v;
}

}  // namespace math
}  // namespace stan
#endif
