#ifndef STAN_MATH_MIX_FUNCTOR_HESSIAN_TIMES_VECTOR_HPP
#define STAN_MATH_MIX_FUNCTOR_HESSIAN_TIMES_VECTOR_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/meta.hpp>

namespace stan {
namespace math {

/**
 * Overload Hessian_times_vector function, under stan/math/mix/functor
 * to handle functions which take in arguments eta, delta, delta_int,
 * and pstream.
 */
template <typename F, typename Eta, require_eigen_t<Eta>* = nullptr>
inline Eigen::VectorXd hessian_times_vector(
    const F& f, const Eigen::VectorXd& x, const Eta& eta,
    const Eigen::VectorXd& delta, const std::vector<int>& delta_int,
    const Eigen::VectorXd& v, std::ostream* pstream = 0) {
  nested_rev_autodiff nested;
  const Eigen::Index x_size = x.size();
  Eigen::Matrix<var, Eigen::Dynamic, 1> x_var = x;
  Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1> x_fvar(x_size);
  for (Eigen::Index i = 0; i < x_size; i++) {
    x_fvar(i) = fvar<var>(x_var(i), v(i));
  }
  fvar<var> fx_fvar = f(x_fvar, eta, delta, delta_int, pstream);
  grad(fx_fvar.d_.vi_);
  return x_var.adj();
}


}  // namespace math
}  // namespace stan

#endif
