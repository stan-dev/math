#ifndef STAN_MATH_MIX_FUNCTOR_HESSIAN_TIMES_VECTOR_HPP
#define STAN_MATH_MIX_FUNCTOR_HESSIAN_TIMES_VECTOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/core.hpp>
#include <stdexcept>
#include <vector>

namespace stan {
namespace math {

template <typename F>
inline void hessian_times_vector(
    const F& f, const Eigen::Matrix<double, Eigen::Dynamic, 1>& x,
    const Eigen::Matrix<double, Eigen::Dynamic, 1>& v, double& fx,
    Eigen::Matrix<double, Eigen::Dynamic, 1>& Hv) {
  using Eigen::Matrix;

  // Run nested autodiff in this scope
  nested_rev_autodiff nested;

  Matrix<var, Eigen::Dynamic, 1> x_var(x.size());
  for (int i = 0; i < x_var.size(); ++i) {
    x_var(i) = x(i);
  }
  var fx_var;
  var grad_fx_var_dot_v;
  gradient_dot_vector(f, x_var, v, fx_var, grad_fx_var_dot_v);
  fx = fx_var.val();
  grad(grad_fx_var_dot_v.vi_);
  Hv.resize(x.size());
  for (int i = 0; i < x.size(); ++i) {
    Hv(i) = x_var(i).adj();
  }
}

template <typename T, typename F, typename EigVec,
          require_eigen_vector_t<EigVec>* = nullptr,
          require_stan_scalar_t<T>* = nullptr>
void hessian_times_vector(const F& f,
                          const Eigen::Matrix<T, Eigen::Dynamic, 1>& x,
                          const EigVec& v, T& fx,
                          Eigen::Matrix<T, Eigen::Dynamic, 1>& Hv) {
  using Eigen::Matrix;
  Matrix<T, Eigen::Dynamic, 1> grad;
  Matrix<T, Eigen::Dynamic, Eigen::Dynamic> H;
  hessian(f, x, fx, grad, H);
  Hv = H * v;
}

/**
 * Overload Hessian_times_vector function, under stan/math/mix/functor
 * to handle functions which take in arguments
 * and pstream.
 */
template <typename F, typename XAdj, typename XVec, typename VVec,
          typename... Args,
          require_all_eigen_vector_t<XAdj, XVec, VVec>* = nullptr>
inline void hessian_times_vector(const F& f, XAdj& x_adj, XVec&& x, VVec&& v,
                                 Args&&... args) {
  nested_rev_autodiff nested;
  const Eigen::Index x_size = x.size();
  Eigen::Matrix<var, Eigen::Dynamic, 1> x_var = std::forward<XVec>(x);
  Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1> x_fvar(x_size);
  for (Eigen::Index i = 0; i < x_size; i++) {
    x_fvar(i) = fvar<var>(x_var(i), v(i));
  }
  fvar<var> fx_fvar = f(x_fvar, args...);
  grad(fx_fvar.d_.vi_);
  x_adj = x_var.adj();
}

}  // namespace math
}  // namespace stan

#endif
