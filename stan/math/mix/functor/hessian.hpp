#ifndef STAN_MATH_MIX_FUNCTOR_HESSIAN_HPP
#define STAN_MATH_MIX_FUNCTOR_HESSIAN_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stdexcept>

namespace stan {
namespace math {

/**
 * Calculate the value, the gradient, and the Hessian,
 * of the specified function at the specified argument in
 * O(N^2) time and O(N^2) space.
 *
 * <p>The functor must implement
 *
 * <code>
 * fvar\<var\>
 * operator()(const
 * Eigen::Matrix\<fvar\<var\>, Eigen::Dynamic, 1\>&)
 * </code>
 *
 * using only operations that are defined for
 * <code>fvar</code> and <code>var</code>.
 *
 * This latter constraint usually
 * requires the functions to be defined in terms of the libraries
 * defined in Stan or in terms of functions with appropriately
 * general namespace imports that eventually depend on functions
 * defined in Stan.
 *
 * @tparam F Type of function
 * @param[in] f Function
 * @param[in] x Argument to function
 * @param[out] fx Function applied to argument
 * @param[out] grad gradient of function at argument
 * @param[out] H Hessian of function at argument
 */
template <typename F>
void hessian(const F& f, const Eigen::Matrix<double, Eigen::Dynamic, 1>& x,
             double& fx, Eigen::Matrix<double, Eigen::Dynamic, 1>& grad,
             Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& H) {
  H.resize(x.size(), x.size());
  grad.resize(x.size());

  // need to compute fx even with size = 0
  if (x.size() == 0) {
    fx = f(x);
    return;
  }
  for (int i = 0; i < x.size(); ++i) {
    // Run nested autodiff in this scope
    nested_rev_autodiff nested;

    Eigen::Matrix<fvar<var>, Eigen::Dynamic, 1> x_fvar(x.size());
    for (int j = 0; j < x.size(); ++j) {
      x_fvar(j) = fvar<var>(x(j), i == j);
    }
    fvar<var> fx_fvar = f(x_fvar);
    grad(i) = fx_fvar.d_.val();
    if (i == 0) {
      fx = fx_fvar.val_.val();
    }
    stan::math::grad(fx_fvar.d_.vi_);
    for (int j = 0; j < x.size(); ++j) {
      H(i, j) = x_fvar(j).val_.adj();
    }
  }
}

}  // namespace math
}  // namespace stan
#endif
