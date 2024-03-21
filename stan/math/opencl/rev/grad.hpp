#ifndef STAN_MATH_OPENCL_REV_FUN_GRAD_HPP
#define STAN_MATH_OPENCL_REV_FUN_GRAD_HPP
#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>

namespace stan {
namespace math {

/**
 * Propagate chain rule to calculate gradients starting from
 * the specified variable.  Resizes the input vector to be the
 * correct size.
 *
 * The grad() function does not itself recover any memory.  use
 * <code>recover_memory()</code> or
 * <code>recover_memory_nested()</code> to recover memory.
 *
 * @param[in] v Value of function being differentiated
 * @param[in] x Variables being differentiated with respect to
 * @param[out] g Gradient, d/dx v, evaluated at x.
 */
inline void grad(var& v, var_value<matrix_cl<double>>& x, Eigen::VectorXd& g) {
  grad(v.vi_);
  g = from_matrix_cl<Eigen::VectorXd>(x.adj());
}

}  // namespace math
}  // namespace stan
#endif
#endif
