#ifndef STAN_MATH_REV_FUN_TO_VAR_HPP
#define STAN_MATH_REV_FUN_TO_VAR_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/core/typedefs.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/prim/fun/typedefs.hpp>
#include <vector>

namespace stan {
namespace math {

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] x A scalar value
 * @return An automatic differentiation variable with the input value.
 */
inline var to_var(double x) { return var(x); }

/**
 * Specialization of to_var for non-const var input
 *
 * @param[in,out] x An automatic differentiation variable.
 * @return The input automatic differentiation variable.
 */
inline var& to_var(var& x) { return x; }

/**
 * Specialization of to_var for const var input
 *
 * @param[in,out] x An automatic differentiation variable.
 * @return The input automatic differentiation variable.
 */
inline const var& to_var(const var& x) { return x; }

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] v A std::vector<double>
 * @return A std::vector<var> with the values set
 */
inline std::vector<var> to_var(const std::vector<double>& v) {
  std::vector<var> var_vector(v.size());
  for (size_t n = 0; n < v.size(); n++) {
    var_vector[n] = v[n];
  }
  return var_vector;
}

/**
 * Specialization of to_var to for const input vector of var
 *
 * Returns a var variable from the input
 *
 * @param[in] v A std::vector<var>
 * @return The input std::vector<var>
 */
inline const std::vector<var>& to_var(const std::vector<var>& v) { return v; }

/**
 * Specialization of to_var to for non-const input vector of var
 *
 * Returns a var variable from the input
 *
 * @param[in] v A std::vector<var>
 * @return The input std::vector<var>
 */
inline std::vector<var>& to_var(std::vector<var>& v) { return v; }

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
 * @param[in,out] m A matrix of automatic differentation variables.
 * @return The input matrix of automatic differentiation variables.
 */
inline matrix_v& to_var(matrix_v& m) { return m; }

/**
 * Specialization of to_var for const matrices of vars
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
 * @param[in,out] v A column vector of automatic differentation variables.
 * @return The input column vector of automatic differentiation variables.
 */
inline const vector_v& to_var(const vector_v& v) { return v; }

/**
 * Specialization of to_var for non-const column vector of vars
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
 * @param[in,out] rv A column vector of automatic differentation variables.
 * @return The input row vector of automatic differentiation variables.
 */
inline const row_vector_v& to_var(const row_vector_v& rv) { return rv; }

/**
 * Specialization of to_var for non-const row vector of vars
 *
 * @param[in,out] rv A column vector of automatic differentation variables.
 * @return The input row vector of automatic differentiation variables.
 */
inline row_vector_v& to_var(row_vector_v& rv) { return rv; }

}  // namespace math
}  // namespace stan
#endif
