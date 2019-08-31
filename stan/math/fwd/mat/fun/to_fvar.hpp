#ifndef STAN_MATH_FWD_MAT_FUN_TO_FVAR_HPP
#define STAN_MATH_FWD_MAT_FUN_TO_FVAR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/scal/fun/to_fvar.hpp>
#include <stan/math/prim/mat/err/check_matching_dims.hpp>

namespace stan {
namespace math {


/**
 * Specialization of to_fvar for matrices of fvars
 *
 * @tparam T The type of the matrix
 * @param[in,out] x A matrix of forward automatic differentation variables.
 * @return The input matrix of forward automatic differentiation variables.
 */
template <typename T, require_eigen_fvar<T>...>
inline auto&& to_fvar(T&& x) {
  return std::forward<T>(x);
}


template <int R, int C>
inline Eigen::Matrix<fvar<double>, R, C> to_fvar(
    const Eigen::Matrix<double, R, C>& x) {
  Eigen::Matrix<fvar<double>, R, C> m_fd(x.rows(), x.cols());
  for (int i = 0; i < x.size(); ++i)
    m_fd(i) = x(i);
  return m_fd;
}

template <typename T, int R, int C>
inline Eigen::Matrix<fvar<T>, R, C> to_fvar(
    const Eigen::Matrix<T, R, C>& val, const Eigen::Matrix<T, R, C>& deriv) {
  check_matching_dims("to_fvar", "value", val, "deriv", deriv);
  Eigen::Matrix<fvar<T>, R, C> ret(val.rows(), val.cols());
  for (int i = 0; i < val.rows(); i++) {
    for (int j = 0; j < val.cols(); j++) {
      ret(i, j).val_ = val(i, j);
      ret(i, j).d_ = deriv(i, j);
    }
  }
  return ret;
}

}  // namespace math
}  // namespace stan
#endif
