#include <stan/math/rev/core.hpp>
#include <stan/math/laplace/lgp_dense_system.hpp>
#include <stan/math/laplace/lgp_dense_newton_solver.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>

// functor for algebraic solver
// struct lgp_functor {
//   template <typename T0, typename T1>
//   inline Eigen::Matrix<typename stan::return_type<T0, T1>::type,
//                        Eigen::Dynamic, 1>
//   operator ()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
//             const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
//             const std::vector<double>& dat,
//             const std::vector<int>& dat_int,
//             std::ostream* pstream__) const {
//     typedef typename stan::return_type<T0, T1>::type scalar;
//     Eigen::Matrix<scalar, Eigen::Dynamic, 1> fgrad;
//     int dim_theta = 2;
//     
//     Eigen::VectorXd n_samples(dim_theta);
//     n_samples(0) = dat[0];
//     n_samples(1) = dat[1];
//     
//     Eigen::VectorXd sums(dim_theta);
//     sums(0) = dat[2];
//     sums(1) = dat[3];
//     
//     return sums - stan::math::elt_multiply(n_samples,
//                                            stan::math::exp(theta))
//       - theta / (phi(0) * phi(0));
//   }
// };

struct lgp_functor {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename stan::return_type<T0, T1>::type,
                       Eigen::Dynamic, 1>
  operator ()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
              const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
              const std::vector<double>& dat,
              const std::vector<int>& dat_int,
              std::ostream* pstream__) const {
    typedef typename stan::return_type<T0, T1>::type scalar;
    using stan::math::lgp_dense_system;
    int dim_theta = 2;
    Eigen::VectorXd n_samples(dim_theta);
    n_samples(0) = dat[0];
    n_samples(1) = dat[1];

    Eigen::VectorXd sums(dim_theta);
    sums(0) = dat[2];
    sums(1) = dat[3];

    lgp_dense_system<scalar> system(phi, n_samples, sums);

    return system.cond_gradient(theta);
  }
};

TEST(laplace, lgp_dense_system) {
  // Since the dense function generalizes the diagonal case (semi-misleadingly
  // labeled conditional), we reproduce the result from the old function on
  // this new case.
  // One caveat: the global parameter in the diagonal case was the standard
  // deviation, but here it is easier to think in terms of variance.
  // Hence the first element of phi is not 2, but 4.
  // R code to generate results are in make_data.R with seed 1954.
  using stan::math::lgp_dense_system;
  using std::cout;
  using std::endl;

  int dim_theta = 2;
  Eigen::VectorXd theta(dim_theta);
  theta << -0.7268203, 1.3347728;

  Eigen::VectorXd phi(2);
  phi << 4, 0;  // global parameters: sigma and rho

  Eigen::VectorXd n_samples(dim_theta);
  for (int i = 0; i < dim_theta; i++) n_samples(i) = 5;
  Eigen::VectorXd sums(dim_theta);
  sums << 3, 10;

  lgp_dense_system<double> system(phi, n_samples, sums);

  // Test evaulation of log densiy
  EXPECT_FLOAT_EQ(6.595955, system.log_density(theta));

  // Test evaluation of the gradient (of the log density)
  Eigen::VectorXd cond_grad = system.cond_gradient(theta);
  EXPECT_FLOAT_EQ(0.7644863, cond_grad(0));
  EXPECT_FLOAT_EQ(-9.3293567, cond_grad(1));

  // Test for the Hessian
  Eigen::MatrixXd cond_hessian = system.cond_hessian(theta);
  EXPECT_FLOAT_EQ(-2.667219, cond_hessian(0, 0));
  EXPECT_FLOAT_EQ(0, cond_hessian(0, 1));
  EXPECT_FLOAT_EQ(0, cond_hessian(1, 0));
  EXPECT_FLOAT_EQ(-19.245664, cond_hessian(1, 1));
}

TEST(laplace, lgp_newton_solver) {
  // Compute mode using Powell's solver and use it as
  // a benchmark.
  using stan::math::var;
  using stan::math::algebra_solver;
  using stan::math::value_of;
  using stan::math::lgp_dense_system;
  using std::cout;
  using std::endl;

  int dim_theta = 2;
  int dim_phi = 2;
  Eigen::VectorXd theta_0(dim_theta);  // initial guess
  theta_0 << -1, 1;
  Eigen::VectorXd n_samples(dim_theta);
  n_samples << 5, 5;
  Eigen::VectorXd sums(dim_theta);
  sums << 3, 10;

  std::vector<double> dat(2 * dim_theta);
  dat[0] = n_samples(0);
  dat[1] = n_samples(1);
  dat[2] = sums(0);
  dat[3] = sums(1);
  std::vector<int> dummy_int;

  // Compute benchmark solution with Powell's method.
  Eigen::MatrixXd solver_gradient(2, 2);
  Eigen::VectorXd powell_solution;

  for (int k = 0; k < dim_theta; k++) {
    var sigma = 1;
    var corr = 0.5;
    Eigen::Matrix<var, Eigen::Dynamic, 1> phi_v(2);
    phi_v << sigma, corr;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
      = algebra_solver(lgp_functor(), theta_0, phi_v, dat, dummy_int);
    powell_solution = value_of(theta);  // a bit redundant...

    AVEC parameters = createAVEC(phi_v(0), phi_v(1));
    VEC g;
    theta(k).grad(parameters, g);
    solver_gradient(k, 0) = g[0];
    solver_gradient(k, 1) = g[1];
  }
  
  // Test newton solver with double arguments.
  Eigen::VectorXd phi_dbl(2);
  phi_dbl << 1, 0.5;
  lgp_dense_system<double> system(phi_dbl, n_samples, sums);
  
  Eigen::VectorXd theta_dbl = lgp_dense_newton_solver(theta_0, system);
  
  EXPECT_FLOAT_EQ(powell_solution(0), theta_dbl(0));
  EXPECT_FLOAT_EQ(powell_solution(1), theta_dbl(1));

}


