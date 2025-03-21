#ifndef STAN_MATH_OPENCL_REV_GRADIENT_HPP
#define STAN_MATH_OPENCL_REV_GRADIENT_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
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
inline void gradient(F&& f, const matrix_cl<double>& x,
              double& fx, matrix_cl<double>& grad_fx) {
  nested_rev_autodiff nested;

  var_value<matrix_cl<double>> x_var(x);
  var fx_var = f(x_var);
  fx = fx_var.val();
  grad(fx_var.vi_);
  grad_fx = x_var.adj();
}


}  // namespace math
}  // namespace stan
#endif
