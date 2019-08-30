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
 * @param[in] m A Matrix with scalars
 * @return A Matrix with automatic differentiation variables
 */
template <typename T, enable_if_eigen<T>..., enable_if_var<scalar_type_decay_t<T>>...>
inline auto&& to_var(T&& m) {
  return std::forward<T>(m);
}

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] m A Matrix with scalars
 * @return A Matrix with automatic differentiation variables
 */
inline matrix_v to_var(const matrix_d& m) {
  matrix_v m_v = m;
  return m_v;
}

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] v A Vector of scalars
 * @return A Vector of automatic differentiation variables with
 *   values of v
 */
inline vector_v to_var(const vector_d& v) {
  vector_v v_v = v;
  return v_v;
}


/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] rv A row vector of scalars
 * @return A row vector of automatic differentation variables with
 *   values of rv.
 */
inline row_vector_v to_var(const row_vector_d& rv) {
  row_vector_v rv_v = rv;
  return rv_v;
}

}  // namespace math
}  // namespace stan
#endif
