#ifndef STAN_MATH_LAPLACE_LG_LOGISTIC_SOLVER_HPP
#define STAN_MATH_LAPLACE_LG_LOGISTIC_SOLVER_HPP

#include <stan/math/rev/mat/functor/kinsol_solve.hpp>

namespace stan {
namespace math {

/**
 * A function on which the Jacobian function can be called to return
 * the derivative of the objective function with respect to the global
 * parameter.
 */
template <typename K>
struct str_phi_sensitivities {
  Eigen::VectorXd theta_;
  std::vector<Eigen::VectorXd> x_;
  K k_f_;

  str_phi_sensitivities(const Eigen::VectorXd& theta,
                        const std::vector<Eigen::VectorXd>& x,
                        const K& k_f) : theta_(theta), x_(x), k_f_(k_f) { }

  template <typename T>
  inline Eigen::Matrix<T, Eigen::Dynamic, 1>
  operator () (const Eigen::Matrix<T, Eigen::Dynamic, 1>& phi) const {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sigma = k_f_(phi, x_);

    return - mdivide_left(Sigma, theta_);
  }
};

/**
 * A functor that returns f, the system function, for
 * a latent Gaussian model with a logistic likelihood. 
 * Note that the number of samples per group, sum per group,
 * (the two being the sufficient stats for theta, the latent
 * gaussian) and the gaussian precision matrix (Q) are 
 * pre-computed and stored in dat.
 */
struct lg_logistic_f {
  inline Eigen::VectorXd
  operator() (const Eigen::VectorXd& theta,
              const Eigen::VectorXd& phi,
              const std::vector<double>& dat,
              const std::vector<int>& dat_int,
              std::ostream* pstream__) const {
    size_t N = theta.size();
    Eigen::Map<const Eigen::VectorXd> n_samples(&dat[0], N);
    Eigen::Map<const Eigen::VectorXd> sums(&dat[N], N);
    Eigen::Map<const Eigen::MatrixXd> Q(&dat[2 * N], N, N);

    return sums - n_samples.cwiseProduct(inv_logit(theta)) - Q * theta;
  }
};

/**
 * A functor that returns the Jacobian of f, the system function,
 * with respect to the unknown theta.
 */
struct lg_logistic_J_f {
  template <typename F>
  inline int
  operator() (const F& f,
              const Eigen::VectorXd& theta,
              const Eigen::VectorXd& phi,
              const std::vector<double>& dat,
              const std::vector<int>& dat_int,
              std::ostream* msgs,
              const double x_sun[], SUNMatrix J) const {
    size_t N = theta.size();

    Eigen::Map<const Eigen::VectorXd> n_samples(&dat[0], N);
    Eigen::Map<const Eigen::MatrixXd> Q(&dat[2 * N], N, N);
    Eigen::Map<const Eigen::VectorXd> x(&x_sun[0], N);

    std::vector<double> jacobian_x = std::vector<double>(N * N);
    Eigen::VectorXd exp_x = exp(x);
    Eigen::VectorXd one = rep_vector(1, N);
    Eigen::VectorXd likelihood_diff
      = elt_divide(exp_x, square(one + exp_x)); 

    Eigen::Map<Eigen::MatrixXd>(&jacobian_x[0], N, N)
      = (- n_samples).cwiseProduct(likelihood_diff).asDiagonal();
    Eigen::Map<Eigen::MatrixXd>(&jacobian_x[0], N, N) -= Q;

    std::move(jacobian_x.begin(), jacobian_x.end(), SM_DATA_D(J));

    return 0;
  }
};

/**
 * Vari class for lg_logistic_solver.
 */
template <typename T>
struct lg_logistic_solver_vari : public vari {
  int phi_size_;
  vari** phi_;
  int theta_size_;
  vari** theta_;
  Eigen::MatrixXd J_;  // CHECK - use double* J_;

  template <typename K>
  lg_logistic_solver_vari(const Eigen::Matrix<T, Eigen::Dynamic, 1>& phi,
                          const std::vector<Eigen::VectorXd>& x,
                          const K& kf,
                          const std::vector<int>& n_samples,
                          const Eigen::VectorXd& theta_dbl)
    : vari(theta_dbl(0)),
      phi_size_(phi.size()),
      phi_(ChainableStack::instance().memalloc_.alloc_array<vari*>(phi_size_)), 
      theta_size_(theta_dbl.size()),
      theta_(ChainableStack::instance().memalloc_.alloc_array<vari*>(
          theta_size_)) {
    for (int i = 0; i < phi_size_; i++) phi_[i] = phi(i).vi_;

    theta_[0] = this;
    for (int i = 0; i < theta_size_; i++)
      theta_[i] = new vari(theta_dbl(i), false);

    // compute derivatives w.r.t phi (using fwd autodiff)
    Eigen::VectorXd phi_dbl = value_of(phi);
    Eigen::MatrixXd phi_sensitivities;
    Eigen::VectorXd dummy;
    str_phi_sensitivities<K> f(theta_dbl, x, kf);
    jacobian_fwd(f, phi_dbl, dummy, phi_sensitivities);

    // compute derivative w.r.t theta (analytically)
    // FIX ME -- do not repeat code. Have one function for Hessian.
    Eigen::VectorXd exp_theta = exp(theta_dbl);
    Eigen::VectorXd one = rep_vector(1, theta_size_);
    Eigen::VectorXd likelihood_diff
      = elt_divide(exp_theta, square(one + exp_theta)); 
    Eigen::MatrixXd theta_sensitivities
      = (- to_vector(n_samples)).cwiseProduct(likelihood_diff).asDiagonal();
    theta_sensitivities -= inverse_spd(kf(phi_dbl, x));

    // apply the implicit function theorem
    J_ = - mdivide_left(theta_sensitivities, phi_sensitivities);
  }

  void chain() {
    for (int j = 0; j < phi_size_; j++)
      for (int i = 0; i < theta_size_; i++)
        phi_[i]->adj_ += theta_[i]->adj_ * J_(i, j);
  }
};

/**
 * Function definition when phi is passed as a vector of double.
 * In future version, the user will specify a functor that returns
 * the covariance matrix.
 *
 * theta_0: intial guess
 * phi: vector of parameters - here contains length scale and std.
 * x: array of covariate vectors (one point per group)
 * k: a functor that returns the latent covariance function.
 * n_samples: number of samples per group
 * sums: sum of responses within one group.
 */
template <typename T, typename K>
Eigen::VectorXd lg_logistic_solver(
  const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta_0,
  const Eigen::VectorXd& phi,
  const std::vector<Eigen::VectorXd>& x,
  const K& k_f,
  const std::vector<int>& n_samples,
  const std::vector<int>& sums,
  double function_tol = 1e-6,
  long int max_num_steps = 100) {  // NOLINT(runtime/int)

  // Construct dat array for system
  int M = theta_0.size();
  std::vector<double> dat (M * (2 + M));
  for (int i = 0; i < M; i++) dat[i] = n_samples[i];
  for (int i = 0; i < M; i++) dat[M + i] = sums[i];
  {
    Eigen::MatrixXd Sigma = k_f(phi, x);
    std::vector<double> Q_array = to_array_1d(inverse_spd(Sigma));
    for (int i = 0; i < M * M; i++) dat[2 * M + i] = Q_array[i];
  }
  std::vector<int> dummy_int;

  // Run Kinsol solver
  Eigen::VectorXd theta_dbl
    = kinsol_solve(lg_logistic_f(), lg_logistic_J_f(), theta_0, phi,
                   dat, dummy_int);
  return theta_dbl;
}

/**
 * Function when phi is passed as a vector of var.
 */
template <typename T1, typename T2, typename K>
Eigen::Matrix<T2, Eigen::Dynamic, 1> lg_logistic_solver(
  const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta_0,
  const Eigen::Matrix<T2, Eigen::Dynamic, 1>& phi,
  const std::vector<Eigen::VectorXd>& x,
  const K& k_f,
  const std::vector<int>& n_samples,
  const std::vector<int>& sums,
  double function_tol = 1e-6,
  long int max_num_steps = 100) {

  Eigen::VectorXd theta_dbl
    = lg_logistic_solver(value_of(theta_0), value_of(phi), x, k_f,
                         n_samples, sums, function_tol, max_num_steps);

  // construct vari
  lg_logistic_solver_vari<T2>* vi0
    = new lg_logistic_solver_vari<T2>(phi, x, k_f, n_samples, theta_dbl);

  Eigen::Matrix<var, Eigen::Dynamic, 1> theta(theta_dbl.size());
  theta(0) = var(vi0->theta_[0]);
  for (int i = 1; i < theta_dbl.size(); i++)
    theta(i) = var(vi0->theta_[i]);

  return theta;
}

}  // namespace math
}  // namespace stan

#endif
