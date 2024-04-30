#ifndef STAN_MATH_MIX_FUNCTOR_HESSIAN_HPP
#define STAN_MATH_MIX_FUNCTOR_HESSIAN_HPP

#include <stan/math/prim/fun/Eigen.hpp>
#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/value_of_rec.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/fun/value_of_rec.hpp>
#include <stdexcept>

namespace stan {
namespace math {

/**
 * Calculate the value, the gradient, and the Hessian,
 * of the specified function at the specified argument in
 * O(N^2) time and O(N^2) space.
 *
 * Instead of returning the full symmetric Hessian, we return the
 * lower-triangular only as a column-major compressed sparse matrix.
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
 * @param[out] H Hessian of function at argument, as a lower-triangular
 *                      compressed sparse matrix
 */
template <typename F>
void hessian(const F& f, const Eigen::VectorXd& x, double& fx,
             Eigen::VectorXd& grad, Eigen::SparseMatrix<double>& H) {
  int d = x.size();
  if (d == 0) {
    fx = value_of_rec(f(x));
    return;
  }

  grad.resize(d);
  H.resize(d, d);
  H.reserve(Eigen::VectorXi::LinSpaced(d, 1, d).reverse());

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
    for (int j = i; j < x.size(); ++j) {
      H.insert(j, i) = x_fvar(j).val_.adj();
    }
  }
  H.makeCompressed();
}

/**
 * Calculate the value, the gradient, and the Hessian,
 * of the specified function at the specified argument in
 * O(N^2) time and O(N^2) space.
 *
 * Overload for returning the Hessian as a symmetric dense matrix.
 *
 * @tparam F Type of function
 * @param[in] f Function
 * @param[in] x Argument to function
 * @param[out] fx Function applied to argument
 * @param[out] grad gradient of function at argument
 * @param[out] H Hessian of function at argument, as a symmetric matrix
 */
template <typename F>
void hessian(const F& f, const Eigen::VectorXd& x, double& fx,
             Eigen::VectorXd& grad, Eigen::MatrixXd& H) {
  Eigen::SparseMatrix<double> hess_sparse;
  hessian(f, x, fx, grad, hess_sparse);

  H = Eigen::MatrixXd(hess_sparse).selfadjointView<Eigen::Lower>();
}

}  // namespace math
}  // namespace stan
#endif
