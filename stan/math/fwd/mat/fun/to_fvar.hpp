#ifndef STAN_MATH_FWD_MAT_FUN_TO_FVAR_HPP
#define STAN_MATH_FWD_MAT_FUN_TO_FVAR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/scal/fun/to_fvar.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/mat/err/check_matching_dims.hpp>

namespace stan {
namespace math {

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] v A Vector of automatic differentiation variables
 * @return A Vector of automatic differentiation variables with
 *   values of v
 */
inline vector_fd& to_fvar(vector_fd& v) { return v; }

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] m A Matrix with automatic differentiation variables.
 * @return A Matrix with automatic differentiation variables.
 */
inline matrix_fd& to_fvar(matrix_fd& m) { return m; }

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] rv A row vector with automatic differentiation variables
 * @return A row vector with automatic differentiation variables
 *    with values of rv.
 */
inline row_vector_fd& to_fvar(row_vector_fd& rv) { return rv; }

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] m A Matrix with scalars
 * @return A Matrix with automatic differentiation variables
 */
inline matrix_fd to_fvar(const matrix_d& m) {
  matrix_fd m_v = m;
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
inline vector_fd to_fvar(const vector_d& v) {
  vector_fd v_v = v;
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
inline row_vector_fd to_fvar(const row_vector_d& rv) {
  row_vector_fd rv_v = rv;
  return rv_v;
}

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] m A Matrix with scalars
 * @return A Matrix with automatic differentiation variables
 */
template <typename T1, int R1, int C1, typename T2, int R2, int C2,
          enable_if_fvar<T2>* = nullptr>
inline Eigen::Matrix<T2, R1, C1> to_fvar(const Eigen::Matrix<T1, R1, C1>& m1,
                                         const Eigen::Matrix<T2, R2, C2>& m2) {
  Eigen::Matrix<T2, R1, C1> m_v = m1;
  return m_v;
}

/**
 * Converts argument to an automatic differentiation variable.
 *
 * Returns a var variable with the input value.
 *
 * @param[in] m A Matrix with scalars
 * @return A Matrix with automatic differentiation variables
 */
template <typename T1, typename T2, int R2, int C2,
          enable_if_fvar<T2>* = nullptr>
inline T1 to_fvar(const T1& m1, const Eigen::Matrix<T2, R2, C2>& m2) {
  T2 m_v(m1);
  return m_v;
}

template <typename T1, int R1, int C1, typename T2,
          typename = std::enable_if_t<ad_promotable<T1, T2>::value>>
inline Eigen::Matrix<T2, R1, C1> to_fvar(const Eigen::Matrix<T1, R1, C1>& m1,
                                         const T2& m2) {
  Eigen::Matrix<T2, R1, C1> m_v = m1;
  return m_v;
}

template <typename T1, typename T2, enable_if_var_or_arithmetic<T2>* = nullptr>
inline T1& to_fvar(T1& m_v, const T2& m_throwaway) {
  return m_v;
}

}  // namespace math
}  // namespace stan
#endif
