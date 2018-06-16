#ifndef STAN_MATH_REV_MAT_FUNCTOR_QUADRATIC_OPTIMIZER_HPP
#define STAN_MATH_REV_MAT_FUNCTOR_QUADRATIC_OPTIMIZER_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/eiquadprog.hpp>
#include <stan/math/rev/mat/functor/jacobian.hpp>
#include <stan/math/rev/mat/fun/dot_product.hpp>
#include <stan/math/prim/mat/fun/dot_product.hpp>
#include <stan/math/prim/mat/fun/value_of.hpp>
#include <stan/math/prim/mat/fun/to_matrix.hpp>
#include <Eigen/dense>
#include <iostream>
#include <string>
#include <vector>

namespace stan {
namespace math {

/**
 * Return analytical form, which will be used to compute Jacobians.
 */
template <typename F1, typename F2, typename F3, typename F4, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
quadratic_optimizer_analytical
  (const F1& fh,
   const F2& fv,
   const F3& fa,
   const F4& fb,
   const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta,
   const std::vector<double>& delta,
   const std::vector<int>& delta_int,
   const Eigen::VectorXd& x,
   double tol = 1e-10) {
  typedef T scalar;

  using Eigen::Matrix;
  using Eigen::Dynamic;

  Matrix<scalar, Dynamic, Dynamic> H = fh(theta, delta, delta_int);
  Matrix<scalar, Dynamic, Dynamic> v = fv(theta, delta, delta_int);
  Matrix<scalar, Dynamic, 1> A = fa(theta, delta, delta_int);
  scalar b = fb(theta, delta, delta_int);
  int size = x.size();

  Matrix<T, Dynamic, 1> x_analytical(x.size());

  Matrix<scalar, Dynamic, 1> sum_x_H(size);
  for (int j = 0; j < size; j++) {
    Matrix<scalar, Dynamic, 1> H_row = H.row(j);
    sum_x_H(j) = dot_product(H_row, x) - x(j) * H(j, j);
  }

  scalar sum_A_H = 0;
  for (int j = 0; j < size; j++) if (x(j) > tol)
    sum_A_H += A(j) * A(j) / H(j, j);

  scalar sum_lambda = 0;
  for (int k = 0; k < size; k++) {
    if (x(k) > tol)
      sum_lambda += A(k) / H(k, k) * (v(k) + 0.5 * sum_x_H(k));
  }
  // Here we write the equality constraint A.x + b = 0,
  // as opposed to A.x = b, hence the substraction.
  sum_lambda -= b;

  Matrix<scalar, Dynamic, 1> x_an(size);
  for (int i = 0; i < x.size(); i++)
    x_an(i) = (1 / H(i, i)) * (A(i) * sum_lambda / sum_A_H
                                 - v(i) - 0.5 * sum_x_H(i));
  return x_an;
}

/**
* A functor which can be passed to the jacobian function.
*/
template <typename Fh, typename Fv, typename Fa, typename Fb>
struct f_theta {
  std::vector<double> delta_;
  std::vector<int> delta_int_;
  Eigen::VectorXd x_;
  Fh fh_;
  Fv fv_;
  Fa fa_;
  Fb fb_;
  double tol_;

  f_theta(const Fh& fh,
          const Fv& fv,
          const Fa& fa,
          const Fb& fb,
          const std::vector<double>& delta,
          const std::vector<int>& delta_int,
          const Eigen::VectorXd& x,
          double tol = 1e-10) : delta_(delta), delta_int_(delta_int),
          x_(x), fh_(fh), fv_(fv), fa_(fa), fb_(fb), tol_(tol) {}

  template<typename T>
  Eigen::Matrix<T, Eigen::Dynamic, 1>
  inline operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta) const {
    return quadratic_optimizer_analytical(fh_, fv_, fa_, fb_,
                                          theta, delta_, delta_int_, x_,
                                          tol_);
  }
};

/**
 * The vari class for the quadratic optimizer.
 */
template <typename F, typename T>
struct quadratic_optimizer_vari : public vari {
  /** vector of parameters */
  vari** theta_;
  /** number of parameters */
  int theta_size_;
  /** vector of solution */
  vari** x_;
  /** number of unknowns */
  int x_size_;
  /** functor for jacobian function */
  F f_;
  /** Jacobian of the solution w.r.t the parameters */
  double* Jx_theta_;

  quadratic_optimizer_vari(const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta,
                           const Eigen::VectorXd& x,
                           const F& f)
    : vari(x(0)),
      theta_(ChainableStack::instance().memalloc_.alloc_array<vari*>(
          theta.size())),
      theta_size_(theta.size()),
      x_(ChainableStack::instance().memalloc_.alloc_array<vari*>(x.size())),
      x_size_(x.size()),
      f_(f),
      Jx_theta_(ChainableStack::instance().memalloc_.alloc_array<double>(
        x_size_ * theta_size_)) {
    using Eigen::Map;
    using Eigen::MatrixXd;

    for (int i = 0; i < theta_size_; i++) theta_[i] = theta(i).vi_;

    x_[0] = this;
    for (int i = 1; i < x_size_; i++) x_[i] = new vari(x(i), false);

    // Compute the Jacobian matrix and store in array.
    Eigen::VectorXd f_eval;  // dummy argument required by jacobian function.
    MatrixXd J;
    jacobian(f_, value_of(theta), f_eval, J);

    Map<MatrixXd>(&Jx_theta_[0], x_size_, theta_size_) = J;
  }

  void chain() {
    for (int j = 0; j < theta_size_; j++)
      for (int i = 0; i < x_size_; i++)
        theta_[i]->adj_ += x_[i]->adj_ * Jx_theta_[j * x_size_ + i];
  }
};

/**
 * Return the solution to the specified quadratic optimization
 * problem. The objective function we wish to optimize has the form
 * 0.5 x'Hx + vx, where x' is the transpose of x, H is a matrix,
 * and v a vector. We minimize this function for x. In addition,
 * the solution is subject to the constraints:
 * (i) ax + b = 0, where a is a vector and b a scalar.
 * (ii) x > 0, that is all solutions are positive.
 * The variables H, v, a, and b are passed as functions of
 * model parameters theta, and fixed variables delta (real data)
 * and delta_int (integer data).
 * 
 * This implementation is based on eiquadprog.
 * See http://www.labri.fr/perso/guenneba/code/QuadProg/
 * 
 * @tparam Fh type of functor which returns H
 * @tparam Fv type of functor which returns v
 * @tparam Fa type of functor which returns a
 * @tparam Fb type of functor which returns b
 * @param[in] fh functor which returns H
 * @param[in] fv functor which returns a
 * @param[in] fa functor which returns v
 * @param[in] fb functor which returns b
 * @param[in] theta vector of parameters for the objective function
 * @param[in] delta continuous data array for objective function
 * @param[in] delta_int integer data array for objective function
 * @param[in] n the number of unknowns we solve for
 * @param[in, out] msgs the print stream for warning messages (??)
 * @param[in] tol the tolerance level; if the optimizer returns a solution
 *                below the this level, this solution is set to 0. This
 *                plays a critical role in the calculation of derivatives.
 */
template <typename Fh, typename Fv, typename Fa, typename Fb>
Eigen::VectorXd
quadratic_optimizer(const Fh& fh, 
                    const Fv& fv, 
                    const Fa& fa, 
                    const Fb& fb,
                    const Eigen::VectorXd& theta,
                    const std::vector<double>& delta,
                    const std::vector<int>& delta_int,
                    int n,
                    std::ostream* msgs = nullptr,
                    double tol = 1e-10) {
  // TODO: make sure msgs stream works and write debuging message.
  using Eigen::VectorXd;

  VectorXd x;
  double f_value;  // declare here to remove warning message

  VectorXd b(1);
  b(0) = fb(theta, delta, delta_int);

  VectorXd A_v = fa(theta, delta, delta_int);
  Eigen::MatrixXd A(1, A_v.size());
  A = to_matrix(A_v);

  f_value = Eigen::solve_quadprog(fh(theta, delta, delta_int),
                                  fv(theta, delta, delta_int), A, b,
                                  Eigen::MatrixXd::Identity(n, n),
                                  Eigen::VectorXd::Zero(n),
                                  x);
  return x;
}

/**
 * Return the solution to the specified quadratic optimization
 * problem. Evaluation of the function with theta a parameter of var.
 */
template <typename Fh, typename Fv, typename Fa, typename Fb, typename T>
Eigen::Matrix<T, Eigen::Dynamic, 1>
quadratic_optimizer(const Fh& fh,
                    const Fv& fv,
                    const Fa& fa,
                    const Fb& fb,
                    const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta,
                    const std::vector<double>& delta,
                    const std::vector<int>& delta_int,
                    int n,
                    std::ostream* msgs = nullptr,  // TODO: print statements
                    double tol = 1e-10) {
  Eigen::VectorXd
    x_dbl = quadratic_optimizer(fh, fv, fa, fb, value_of(theta), 
                                delta, delta_int, n, msgs);

  // Construct vari
  typedef f_theta<Fh, Fv, Fa, Fb> f_vari;
  quadratic_optimizer_vari<f_vari, T>* vi0
    = new quadratic_optimizer_vari<f_vari, T>(theta, x_dbl,
      f_vari(fh, fv, fa, fb, delta, delta_int, x_dbl, tol));

  Eigen::Matrix<var, Eigen::Dynamic, 1> x(x_dbl.size());
  x(0) = var(vi0->x_[0]);
  for (int i = 1; i < x.size(); i++)
    x(i) = var(vi0->x_[i]);

  return x;
}

}  // namespace math
}  // namespace stan

#endif
