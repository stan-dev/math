#ifndef STAN_MATH_REV_FUNCTOR_GRADIENT_HPP
#define STAN_MATH_REV_FUNCTOR_GRADIENT_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/Eigen.hpp>
#include <stdexcept>
#include <vector>

namespace stan {
namespace math {

/**
 * Calculate the value and the gradient of the specified function
 * at the specified argument.
 *
 * <p>The functor must implement
 *
 * <code>
 * var
 * operator()(const
 * Eigen::Matrix<var, Eigen::Dynamic, 1>&)
 * </code>
 *
 * using only operations that are defined for
 * <code>var</code>.  This latter constraint usually
 * requires the functions to be defined in terms of the libraries
 * defined in Stan or in terms of functions with appropriately
 * general namespace imports that eventually depend on functions
 * defined in Stan.
 *
 * The evaluated gradient is stored into a
 * <code>Eigen::VectorXd</code> named <code>grad_fx</code>.
 *
 * <p>Time and memory usage is on the order of the size of the
 * fully unfolded expression for the function applied to the
 * argument, independently of dimension.
 *
 * @tparam F Type of function
 * @param[in] f Function
 * @param[in] x Argument to function
 * @param[out] fx Function applied to argument
 * @param[out] grad_fx Gradient of function at argument
 */
template <typename F>
void gradient(const F& f, const Eigen::Matrix<double, Eigen::Dynamic, 1>& x,
              double& fx, Eigen::Matrix<double, Eigen::Dynamic, 1>& grad_fx) {
  nested_rev_autodiff nested;

  Eigen::Matrix<var, Eigen::Dynamic, 1> x_var(x);
  var fx_var = f(x_var);
  fx = fx_var.val();
  grad_fx.resize(x.size());
  grad(fx_var.vi_);
  grad_fx = x_var.adj();
}

/**
 * Calculate the value and the gradient of the specified function
 * at the specified argument.
 *
 * <p>The functor must implement
 *
 * <code>
 * var
 * operator()(const
 * Eigen::Matrix<var, Eigen::Dynamic, 1>&)
 * </code>
 *
 * using only operations that are defined for
 * <code>var</code>.  This latter constraint usually
 * requires the functions to be defined in terms of the libraries
 * defined in Stan or in terms of functions with appropriately
 * general namespace imports that eventually depend on functions
 * defined in Stan.
 *
 * The evaluated gradient is stored into the object whose data
 * begins at <code>*first_grad_fx</code> and ends at
 * <code>*last_grad_fx</code>.  The caller is responsible for
 * ensuring the size of the object pointed to by
 * <code>first_grad_fx</code> matches the size of the argument
 * <code>x</code>.
 *
 * <p>Time and memory usage is on the order of the size of the
 * fully unfolded expression for the function applied to the
 * argument, independently of dimension.
 *
 * @tparam F Type of function
 * @param[in] f Function
 * @param[in] x Argument to function
 * @param[out] fx Function applied to argument
 * @param[out] grad_fx Gradient of function at argument
 */
template <typename F, typename InputIt>
void gradient(const F& f, const Eigen::Matrix<double, Eigen::Dynamic, 1>& x,
              double& fx, InputIt first_grad_fx, InputIt last_grad_fx) {
  nested_rev_autodiff nested;

  Eigen::Matrix<var, Eigen::Dynamic, 1> x_var(x);
  var fx_var = f(x_var);
  fx = fx_var.val();
  grad(fx_var.vi_);
  for (size_t i = 0;
       first_grad_fx != last_grad_fx && i < x.size();
       ++first_grad_fx, ++i) {
    *first_grad_fx = x_var(i).adj();
  }
}

}  // namespace math
}  // namespace stan
#endif
