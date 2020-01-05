#ifndef STAN_MATH_FWD_FUN_TO_FVAR_HPP
#define STAN_MATH_FWD_FUN_TO_FVAR_HPP

#include <stan/math/prim/err.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T>
inline fvar<T> to_fvar(const T& x) {
  return fvar<T>(x);
}

/**
 * Specialization of to_fvar for const fvars
 *
 * @param[in,out] x A forward automatic differentation variables.
 * @return The input forward automatic differentiation variables.
 */
template <typename T>
inline const fvar<T>& to_fvar(const fvar<T>& x) {
  return x;
}

/**
 * Specialization of to_fvar for non-const fvars
 *
 * @param[in,out] x A forward automatic differentation variables.
 * @return The input forward automatic differentiation variables.
 */
template <typename T>
inline fvar<T>& to_fvar(fvar<T>& x) {
  return x;
}

template <typename T>
inline std::vector<fvar<T>> to_fvar(const std::vector<T>& v) {
  std::vector<fvar<T>> x(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    x[i] = T(v[i]);
  }
  return x;
}

template <typename T>
inline std::vector<fvar<T>> to_fvar(const std::vector<T>& v,
                                    const std::vector<T>& d) {
  std::vector<fvar<T>> x(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    x[i] = fvar<T>(v[i], d[i]);
  }
  return x;
}

/**
 * Specialization of to_fvar for const fvar input
 *
 * @tparam The inner type of the fvar.
 * @param[in,out] v A vector of forward automatic differentiation variable.
 * @return The input vector of forward automatic differentiation variable.
 */
template <typename T>
inline const std::vector<fvar<T>>& to_fvar(const std::vector<fvar<T>>& v) {
  return v;
}

/**
 * Specialization of to_fvar for non-const fvar input
 *
 * @tparam The inner type of the fvar.
 * @param[in,out] v A vector of forward automatic differentiation variable.
 * @return The input vector of forward automatic differentiation variable.
 */
template <typename T>
inline std::vector<fvar<T>>& to_fvar(std::vector<fvar<T>>& v) {
  return v;
}

/**
 * Specialization of to_fvar for const matrices of fvars
 *
 * @param[in,out] m A matrix of forward automatic differentation variables.
 * @return The input matrix of forward automatic differentiation variables.
 */
template <int R, int C, typename T>
inline const Eigen::Matrix<T, R, C>& to_fvar(const Eigen::Matrix<T, R, C>& m) {
  return m;
}

/**
 * Specialization of to_fvar for non-const matrices of fvars
 *
 * @param[in,out] m A matrix of forward automatic differentation variables.
 * @return The input matrix of forward automatic differentiation variables.
 */
template <int R, int C, typename T>
inline Eigen::Matrix<T, R, C>& to_fvar(Eigen::Matrix<T, R, C>& m) {
  return m;
}

template <int R, int C>
inline Eigen::Matrix<fvar<double>, R, C> to_fvar(
    const Eigen::Matrix<double, R, C>& m) {
  Eigen::Matrix<fvar<double>, R, C> m_fd(m.rows(), m.cols());
  for (int i = 0; i < m.size(); ++i) {
    m_fd(i) = m(i);
  }
  return m_fd;
}

template <typename T, int R, int C>
inline Eigen::Matrix<fvar<T>, R, C> to_fvar(
    const Eigen::Matrix<T, R, C>& val, const Eigen::Matrix<T, R, C>& deriv) {
  check_matching_dims("to_fvar", "value", val, "deriv", deriv);
  Eigen::Matrix<fvar<T>, R, C> ret(val.rows(), val.cols());
  for (int j = 0; j < val.cols(); j++) {
    for (int i = 0; i < val.rows(); i++) {
      ret(i, j).val_ = val(i, j);
      ret(i, j).d_ = deriv(i, j);
    }
  }
  return ret;
}

}  // namespace math
}  // namespace stan
#endif
