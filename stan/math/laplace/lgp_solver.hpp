#ifndef STAN_MATH_LAPLACE_LGP_SOLVER_HPP
#define STAN_MATH_LAPLACE_LGP_SOLVER_HPP

#include <stan/math/rev/mat/functor/kinsol_solve.hpp>
#include <stan/math/laplace/lgp_dense_system.hpp>
#include <stan/math/laplace/lgp_dense_newton_solver.hpp>

namespace stan {
namespace math {

/**
 * A functor that returns f, the system function, for
 * a latent Gaussian Poisson process. Note that the
 * number of samples per group, sum per group,
 * and Gaussian precision matrix (Q) are pre-computed.
 * 
 * CHECK - should we sums and n_samples as vector<int>? 
 */
struct lgp_f {
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

    return sums - n_samples.cwiseProduct(exp(theta)) - Q * theta;
  }
};

/**
 * A functor that returns the Jacobian of f, the system function,
 * with respect to the unknown theta.
 * CHECK -- do I need that strict a signature, or can I relax it?
 */
struct lgp_J_f {
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
    Eigen::Map<Eigen::MatrixXd>(&jacobian_x[0], N, N)
      = (- n_samples).cwiseProduct(exp(x)).asDiagonal();
    Eigen::Map<Eigen::MatrixXd>(&jacobian_x[0], N, N) -= Q;

    std::move(jacobian_x.begin(), jacobian_x.end(), SM_DATA_D(J));

    return 0;
  }
};

/**
 * Function definition when phi is passed as a vector of double.
 */
template <typename T>
Eigen::VectorXd lgp_solver(
  const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta_0,
  const Eigen::VectorXd& phi,
  const std::vector<int>& n_samples,
  const std::vector<int>& sums,
  double function_tol = 1e-6,
  long int max_num_steps = 100) {  // NOLINT(runtime/int)

  bool space_matters = 1;
  lgp_dense_system<double> system(phi,
                                  to_vector(n_samples),
                                  to_vector(sums),
                                  space_matters);

  std::vector<int> dummy_int;
  Eigen::VectorXd theta_dbl
    = kinsol_solve(lgp_f(), lgp_J_f(), theta_0, phi,
                   system.get_dat(), dummy_int);
  return theta_dbl;
}

/**
 * Overload function to directly take in system as an argument.
 */
template <typename T>
Eigen::VectorXd lgp_solver(
    const Eigen::Matrix<T, Eigen::Dynamic, 1>& theta_0,
    const Eigen::VectorXd& phi,
    const lgp_dense_system<double>& system,
    double tol = 1e-3,
    long int max_num_steps = 100)  {
  std::vector<int> dummy_int;
  Eigen::VectorXd theta_dbl
    = kinsol_solve(lgp_f(), lgp_J_f(), theta_0, phi,
                   system.get_dat(), dummy_int);
  return theta_dbl;
}

/**
 * Function when phi is passed as a vector of var.
 */
template <typename T1, typename T2>
Eigen::Matrix<T2, Eigen::Dynamic, 1> lgp_solver(
  const Eigen::Matrix<T1, Eigen::Dynamic, 1>& theta_0,
  const Eigen::Matrix<T2, Eigen::Dynamic, 1>& phi,
  const std::vector<int>& n_samples,
  const std::vector<int>& sums,
  double function_tol = 1e-6,
  long int max_num_steps = 100) {

  bool space_matters = true;
  lgp_dense_system<double> system(value_of(phi),
                                  to_vector(n_samples),
                                  to_vector(sums),
                                  space_matters);

  Eigen::VectorXd theta_dbl
    = lgp_solver(value_of(theta_0), value_of(phi), system,
                 function_tol, max_num_steps);

  // construct vari
  lgp_dense_newton_solver_vari<T2>* vi0
    = new lgp_dense_newton_solver_vari<T2>(phi, system, theta_dbl);

  Eigen::Matrix<var, Eigen::Dynamic, 1> theta(theta_dbl.size());
  theta(0) = var(vi0->theta_[0]);
  for (int i = 1; i < theta_dbl.size(); i++)
    theta(i) = var(vi0->theta_[i]);

  return theta;
}

}  // namespace math
}  // namespace stan

#endif
