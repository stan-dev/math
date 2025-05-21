#ifndef STAN_MATH_FWD_FUN_TO_FVAR_HPP
#define STAN_MATH_FWD_FUN_TO_FVAR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/meta.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/err.hpp>
#include <vector>

namespace stan {
namespace math {

template <typename T, require_stan_scalar_t<T>* = nullptr,
          require_not_fvar_t<T>* = nullptr>
inline fvar<T> to_fvar(const T& x) {
  return fvar<T>(x);
}

/**
 * Specialization of to_fvar for [containers of] fvars
 *
 * @param[in,out] x A forward automatic differentation variables.
 * @return The input forward automatic differentiation variables.
 */
template <typename T, require_fvar_t<scalar_type_t<T>>* = nullptr>
inline T&& to_fvar(T&& x) {
  return std::forward<T>(x);
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

template <typename T, require_eigen_t<T>* = nullptr,
          require_not_eigen_vt<is_fvar, T>* = nullptr>
inline promote_scalar_t<fvar<value_type_t<T>>, T> to_fvar(const T& m) {
  promote_scalar_t<fvar<value_type_t<T>>, T> m_fd(m.rows(), m.cols());
  m_fd.val() = m;
  m_fd.d() = plain_type_t<T>::Constant(m.rows(), m.cols(), 0);
  return m_fd;
}

template <typename T1, typename T2, require_all_eigen_t<T1, T2>* = nullptr,
          require_vt_same<T1, T2>* = nullptr>
inline promote_scalar_t<fvar<value_type_t<T1>>, T1> to_fvar(const T1& val,
                                                            const T2& deriv) {
  check_matching_dims("to_fvar", "value", val, "deriv", deriv);
  promote_scalar_t<fvar<value_type_t<T1>>, T1> ret(val.rows(), val.cols());
  ret.val() = val;
  ret.d() = deriv;
  return ret;
}

}  // namespace math
}  // namespace stan
#endif
