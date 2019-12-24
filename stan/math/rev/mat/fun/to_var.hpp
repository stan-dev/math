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
inline matrix_v to_var(const matrix_d& m) {
  matrix_v m_v = m;
  return m_v;
}

/**
 * Specialization of to_var for non-const matrices of vars
 *
 *
 * @param[in,out] m A matrix of automatic differentation variables.
 * @return The input matrix of automatic differentiation variables.
 */
inline matrix_v& to_var(matrix_v& m) { return m; }

/**
 * Specialization of to_var for const matrices of vars
 *
 *
 * @param[in,out] m A matrix of automatic differentation variables.
 * @return The input matrix of automatic differentiation variables.
 */
inline const matrix_v& to_var(const matrix_v& m) { return m; }

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
 * Specialization of to_var for const column vector of vars
 *
 *
 * @param[in,out] v A column vector of automatic differentation variables.
 * @return The input column vector of automatic differentiation variables.
 */
inline const vector_v& to_var(const vector_v& v) { return v; }

/**
 * Specialization of to_var for non-const column vector of vars
 *
 *
 * @param[in,out] v A column vector of automatic differentation variables.
 * @return The input column vector of automatic differentiation variables.
 */
inline vector_v& to_var(vector_v& v) { return v; }

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

/**
 * Specialization of to_var for const row vector of vars
 *
 *
 * @param[in,out] rv A column vector of automatic differentation variables.
 * @return The input row vector of automatic differentiation variables.
 */
inline const row_vector_v& to_var(const row_vector_v& rv) { return rv; }

/**
 * Specialization of to_var for non-const row vector of vars
 *
 *
 * @param[in,out] rv A column vector of automatic differentation variables.
 * @return The input row vector of automatic differentiation variables.
 */
inline row_vector_v& to_var(row_vector_v& rv) { return rv; }

}  // namespace math
}  // namespace stan
#endif
