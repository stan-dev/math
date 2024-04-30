#ifndef STAN_MATH_FWD_FUNCTOR_HESSIAN_HPP
#define STAN_MATH_FWD_FUNCTOR_HESSIAN_HPP

#include <stan/math/fwd/core.hpp>
#include <stan/math/fwd/fun/value_of.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Calculate the value, the gradient, and the Hessian,
 * of the specified function at the specified argument in
 * time O(N^3) time and O(N^2) space.  The advantage over the
 * mixed definition, which is faster for Hessians, is that this
 * version is itself differentiable.
 *
 * Instead of returning the full symmetric Hessian, we return the
 * lower-triangular only as a column-major compressed sparse matrix.
 *
 * <p>The functor must implement
 *
 * <code>
 * fvar\<fvar\<T\> \>
 * operator()(const
 * Eigen::Matrix\<fvar\<fvar\<T\> \>, Eigen::Dynamic, 1\>&)
 * </code>
 *
 * using only operations that are defined for the argument type.
 *
 * This latter constraint usually requires the functions to be
 * defined in terms of the libraries defined in Stan or in terms
 * of functions with appropriately general namespace imports that
 * eventually depend on functions defined in Stan.
 *
 * @tparam T type of elements in the vector and matrix
 * @tparam F type of function
 * @param[in] f Function
 * @param[in] x Argument to function
 * @param[out] fx Function applied to argument
 * @param[out] grad gradient of function at argument
 * @param[out] H Hessian of function at argument, as a lower-triangular
 *                      compressed sparse matrix
 */
template <typename T, typename F>
void hessian(const F& f, const Eigen::Matrix<T, Eigen::Dynamic, 1>& x, T& fx,
             Eigen::Matrix<T, Eigen::Dynamic, 1>& grad,
             Eigen::SparseMatrix<T>& H) {
  int d = x.size();
  if (d == 0) {
    fx = value_of(f(x));
    return;
  }

  H.resize(d, d);
  H.reserve(Eigen::VectorXi::LinSpaced(d, 1, d).reverse());
  grad.resize(d);

  Eigen::Matrix<fvar<fvar<T> >, Eigen::Dynamic, 1> x_fvar(d);
  for (int i = 0; i < d; ++i) {
    for (int j = i; j < d; ++j) {
      for (int k = 0; k < d; ++k) {
        x_fvar(k) = fvar<fvar<T> >(fvar<T>(x(k), j == k), fvar<T>(i == k, 0));
      }
      fvar<fvar<T> > fx_fvar = f(x_fvar);
      if (j == 0) {
        fx = fx_fvar.val_.val_;
      }
      if (i == j) {
        grad(i) = fx_fvar.d_.val_;
      }
      H.insert(j, i) = fx_fvar.d_.d_;
    }
  }
  H.makeCompressed();
}

/**
 * Calculate the value, the gradient, and the Hessian,
 * of the specified function at the specified argument in
 * time O(N^3) time and O(N^2) space.  The advantage over the
 * mixed definition, which is faster for Hessians, is that this
 * version is itself differentiable.
 *
 * Overload for returning the Hessian as a symmetric dense matrix.
 *
 * @tparam T type of elements in the vector and matrix
 * @tparam F type of function
 * @param[in] f Function
 * @param[in] x Argument to function
 * @param[out] fx Function applied to argument
 * @param[out] grad gradient of function at argument
 * @param[out] H Hessian of function at argument, as a symmetric matrix
 */
template <typename T, typename F>
void hessian(const F& f, const Eigen::Matrix<T, Eigen::Dynamic, 1>& x, T& fx,
             Eigen::Matrix<T, Eigen::Dynamic, 1>& grad,
             Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& H) {
  Eigen::SparseMatrix<T> hess_sparse;
  hessian(f, x, fx, grad, hess_sparse);

  H = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>(hess_sparse)
          .template selfadjointView<Eigen::Lower>();
}

}  // namespace math
}  // namespace stan
#endif
