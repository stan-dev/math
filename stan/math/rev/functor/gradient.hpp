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
 * @tparam EigVec Type of Eigen vector
 * @tparam InputIt must meet the requirements of
 * [LegacyInputIterator](https://en.cppreference.com/w/cpp/named_req/InputIterator).
 * @param[in] f Function
 * @param[in] x Argument to function
 * @param[out] fx Function applied to argument
 * @param[out] first_grad_fx First element of gradient of function at argument
 * @param[out] last_grad_fx Last element of gradient of function at argument
 * @throw std::invalid_argument if the iterator isn't the right size
 * to hold the gradients
 */
template <typename F, typename EigVec, typename InputIt,
          require_eigen_vector_vt<std::is_arithmetic, EigVec>* = nullptr>
void gradient(const F& f, const EigVec& x, double& fx, InputIt first_grad_fx,
              InputIt last_grad_fx) {
  nested_rev_autodiff nested;

  if (last_grad_fx - first_grad_fx != x.size()) {
    std::stringstream s;
    s << "gradient(): iterator and gradient different sizes; iterator size = "
      << last_grad_fx - first_grad_fx << "; grad size = " << x.size()
      << std::endl;
    throw std::invalid_argument(s.str());
  }

  Eigen::Matrix<var, Eigen::Dynamic, 1> x_var(x);
  var fx_var = f(x_var);
  fx = fx_var.val();
  grad(fx_var.vi_);
  for (Eigen::VectorXd::Index i = 0; i < x_var.size(); ++i) {
    *first_grad_fx++ = x_var.coeff(i).adj();
  }
}

}  // namespace math
}  // namespace stan
#endif
