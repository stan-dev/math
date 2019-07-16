#include <stan/math.hpp>
#include <stan/math/laplace/lg_logistic_solver.hpp>
#include <stan/math/laplace/lgp_density.hpp>
#include <stan/math/laplace/laplace_marginal.hpp>
#include <test/unit/math/laplace/lgp_utility.hpp>

#include <test/unit/math/rev/mat/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

TEST(laplace, lg_logistic_solver) {
  using stan::math::lg_logistic_solver;
  using stan::math::lg_logistic_f;
  using stan::math::to_array_1d;
  using stan::math::inverse_spd;
  using stan::math::var;
  using stan::math::diff_logistic_log;
  using stan::math::to_vector;

  // global parameter: variance and length scale
  int dim_phi = 2;
  Eigen::VectorXd phi(dim_phi);
  phi << 1, 1;

  int dim_theta = 2;
  Eigen::VectorXd theta_0(dim_theta);
  theta_0 << 0, 0;

  int dim_x = 2;
  std::vector<Eigen::VectorXd> x(dim_theta);
  Eigen::VectorXd x_0(2);
  x_0 <<  0.05100797, 0.16086164;
  Eigen::VectorXd x_1(2);
  x_1 << -0.59823393, 0.98701425;
  x[0] = x_0;
  x[1] = x_1;

  std::vector<int> n_samples = {1, 1};
  std::vector<int> sums = {0, 1};

  // Check Kernel function
  squared_kernel_functor K;
  // std::cout << "kernel function: " << std::endl
  //           << K(phi, x) << std::endl;

  /////////////////////////////////////////////////////////////////////////////
  // Test lg solver.

  // Run solver
  Eigen::VectorXd
    theta = lg_logistic_solver(theta_0, phi, x, squared_kernel_functor(),
                               n_samples, sums);

  // Evaluate target function
  int M = theta_0.size();
  std::vector<double> dat (M * (2 + M));
  for (int i = 0; i < M; i++) dat[i] = n_samples[i];
  for (int i = 0; i < M; i++) dat[M + i] = sums[i];
  Eigen::MatrixXd Sigma = K(phi, x);
  std::vector<double> Q_array = to_array_1d(inverse_spd(Sigma));
  for (int i = 0; i < M * M; i++) dat[2 * M + i] = Q_array[i];
  std::vector<int> dummy_int;
  lg_logistic_f f;
  Eigen::VectorXd evaluation = f(theta, phi, dat, dummy_int, 0);

  EXPECT_NEAR(0, evaluation(0), 1e-5);
  EXPECT_NEAR(0, evaluation(1), 1e-5);

  // Compute the marginal density and its derivatives
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi_v = phi;
  Eigen::Matrix<var, Eigen::Dynamic, 1>
    theta_v = lg_logistic_solver(theta_0, phi_v, x, squared_kernel_functor(),
                                 n_samples, sums);


  diff_logistic_log diff_likelihood(to_vector(n_samples), to_vector(sums));

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
    Sigma_v = K(phi_v, x, 0);
  var target = diff_likelihood.log_likelihood(theta_v);
  Eigen::Matrix<var, Eigen::Dynamic, 1> gradient, hessian;
  diff_likelihood.diff(theta_v, gradient, hessian);

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
    first_matrix = Eigen::MatrixXd::Identity(dim_theta, dim_theta)
    - diag_post_multiply(Sigma_v, hessian);

  target += - 0.5 * (log_determinant(first_matrix) 
                       + dot_product(theta_v, mdivide_left(Sigma_v, theta)));
  VEC g;
  AVEC parm_vec = createAVEC(phi_v(0), phi_v(1));
  target.grad(parm_vec, g);

  std::cout << "mode: " << theta.transpose() << std::endl 
            << "target: " << target.val() << std::endl
            << "lg grad: " << g[0] << " " << g[1] << std::endl;
            // << "Sigma: " << value_of(Sigma_v) << std::endl;

  /////////////////////////////////////////////////////////////////////////////
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi_v2 = phi;

  var target_laplace = laplace_marginal_density(theta, phi_v2, x,
    diff_logistic_log(to_vector(n_samples), to_vector(sums)),
    squared_kernel_functor(), 1e-3, 100);

  VEC g2;
  AVEC parm_vec2 = createAVEC(phi_v2(0), phi_v2(1));
  target_laplace.grad(parm_vec2, g2);

  std::cout << "target: " << target_laplace.val() << std::endl
            << "laplace grad: " << g2[0] << " " << g[1] << std::endl;
  
  // Use finite diff
  double diff = 1e-7;
  Eigen::VectorXd phi_1l = phi;
  Eigen::VectorXd phi_1u = phi;
  phi_1l(0) = phi(0) - diff;
  phi_1u(0) = phi(0) + diff;
  
  double
  target_1u = laplace_marginal_density(theta, phi_1u, x,
                 diff_logistic_log(to_vector(n_samples), to_vector(sums)),
                 squared_kernel_functor(), 1e-3, 100);

  double
  target_1l = laplace_marginal_density(theta, phi_1l, x,
                 diff_logistic_log(to_vector(n_samples), to_vector(sums)),
                 squared_kernel_functor(), 1e-3, 100);

  Eigen::VectorXd phi_2l = phi;
  Eigen::VectorXd phi_2u = phi;
  phi_2l(1) = phi(1) - diff;
  phi_2u(1) = phi(1) + diff;

  double
  target_2u = laplace_marginal_density(theta, phi_2u, x,
                 diff_logistic_log(to_vector(n_samples), to_vector(sums)),
                 squared_kernel_functor(), 1e-3, 100);

  double
  target_2l = laplace_marginal_density(theta, phi_2l, x,
                 diff_logistic_log(to_vector(n_samples), to_vector(sums)),
                 squared_kernel_functor(), 1e-3, 100);

  VEC g3(2);
  g3[0] = (target_1u - target_1l) / (2 * diff);
  g3[1] = (target_2u - target_2l) / (2 * diff);
  std::cout << "finite diff: " << g3[0] << " " << g3[1] << std::endl;
}
    
    
    

TEST(laplace, likelihood_differentiation) {
  using stan::math::diff_logistic_log;
  using stan::math::var;

  double test_tolerance = 2e-4;

  Eigen::VectorXd theta(2);
  theta << -2.45809, -3.6127;
  Eigen::VectorXd y(2), n_samples(2);
  y << 1, 0;
  n_samples << 1, 1;
  Eigen::Matrix<var, Eigen::Dynamic, 1> theta_v = theta;

  diff_logistic_log diff_functor(n_samples, y);
  double log_density = diff_functor.log_likelihood(theta);
  Eigen::VectorXd gradient, hessian;
  diff_functor.diff(theta, gradient, hessian);
  Eigen::VectorXd third_tensor = diff_functor.third_diff(theta);

  EXPECT_NEAR(-2.566843, log_density, test_tolerance);

  // finite diff calculations for first-order derivatives
  double diff = 1e-12;
  Eigen::VectorXd theta_1u = theta;
  Eigen::VectorXd theta_1l = theta;
  Eigen::VectorXd theta_2u = theta;
  Eigen::VectorXd theta_2l = theta;
  theta_1u(0) = theta(0) + diff;
  theta_1l(0) = theta(0) - diff;
  theta_2u(1) = theta(1) + diff;
  theta_2l(1) = theta(1) - diff;
  double diff_1 = (diff_functor.log_likelihood(theta_1u)
                     - diff_functor.log_likelihood(theta_1l)) / (2 * diff);
  double diff_2 = (diff_functor.log_likelihood(theta_2u)
                     - diff_functor.log_likelihood(theta_2l)) / (2 * diff);

  EXPECT_NEAR(diff_1, gradient(0), test_tolerance);
  EXPECT_NEAR(diff_2, gradient(1), test_tolerance);

  // finite diff calculation for second-order derivatives
  Eigen::VectorXd gradient_1u, gradient_1l, hessian_1u, hessian_1l,
    gradient_2u, gradient_2l, hessian_2u, hessian_2l;
  diff_functor.diff(theta_1u, gradient_1u, hessian_1u);
  diff_functor.diff(theta_1l, gradient_1l, hessian_1l);
  diff_functor.diff(theta_2u, gradient_2u, hessian_2u);
  diff_functor.diff(theta_2l, gradient_2l, hessian_2l);

  double diff_grad_1 = (gradient_1u(0) - gradient_1l(0)) / (2 * diff);
  double diff_grad_2 = (gradient_2u(1) - gradient_2l(1)) / (2 * diff);
  
  EXPECT_NEAR(diff_grad_1, hessian(0), test_tolerance);
  EXPECT_NEAR(diff_grad_2, hessian(1), test_tolerance);
  
  // finite diff calculation for third-order derivatives
  double diff_hess_1 = (hessian_1u(0) - hessian_1l(0)) / (2 * diff);
  double diff_hess_2 = (hessian_2u(1) - hessian_2l(1)) / (2 * diff);

  EXPECT_NEAR(diff_hess_1, third_tensor(0), test_tolerance);
  EXPECT_NEAR(diff_hess_2, third_tensor(1), test_tolerance);
}



TEST(laplace, lg_logistic_solver2) {
  using stan::math::var;
  using stan::math::lg_logistic_solver;
  using stan::math::lg_logistic_f;
  using stan::math::diff_logistic_log;
  using stan::math::to_vector;
  using stan::math::to_array_1d;
  using stan::math::inverse_spd;
  using stan::math::log_determinant;
  using stan::math::mdivide_left;
  using stan::math::start_nested;
  using stan::math::recover_memory_nested;

  int dim_theta = 500;
  int n_observations = 500;
  std::string data_directory = "test/unit/math/laplace/aki_synth_data/";
  std::vector<double> x1(dim_theta), x2(dim_theta);
  std::vector<int> y(n_observations);
  read_in_data(dim_theta, n_observations, data_directory, x1, x2, y);

  // Look a some of the data.
  // std::cout << "x_1: " << x1[0] << " " << x2[0] << std::endl
  //           << "x_2: " << x1[1] << " " << x2[1] << std::endl
  //           << "y_1: " << y[0] << " y_2: " << y[1] << std::endl; 

  int dim_x = 2;
  std::vector<Eigen::VectorXd> x(dim_theta);
  for (int i = 0; i < dim_theta; i++) {
    Eigen::VectorXd coordinate(dim_x);
    coordinate << x1[i], x2[i];
    x[i] = coordinate;
  }

  Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);
  Eigen::VectorXd phi(2);
  phi << 1.6, 1;  // standard deviation, length scale
  std::vector<int> n_samples = stan::math::rep_array(1, dim_theta);

  auto start_optimization = std::chrono::system_clock::now();
  Eigen::VectorXd
    theta = lg_logistic_solver(theta_0, phi, x, squared_kernel_functor(),
                               n_samples, y, 1, 100);
  auto end_optimization = std::chrono::system_clock::now();
  std::chrono::duration<double>
    elapsed_time_optimization = end_optimization - start_optimization;

  // Evaluate the norm of the target function
  int M = theta_0.size();
  std::vector<double> dat (M * (2 + M));
  for (int i = 0; i < M; i++) dat[i] = n_samples[i];
  for (int i = 0; i < M; i++) dat[M + i] = y[i];
  squared_kernel_functor K;
  Eigen::MatrixXd Sigma = K(phi, x);
  std::vector<double> Q_array = to_array_1d(inverse_spd(Sigma));
  for (int i = 0; i < M * M; i++) dat[2 * M + i] = Q_array[i];
  std::vector<int> dummy_int;
  lg_logistic_f f;
  Eigen::VectorXd evaluation = f(theta, phi, dat, dummy_int, 0);

  std::cout << "target norm: " << evaluation.norm() << std::endl
            << "solving time: " << elapsed_time_optimization.count()
            << std::endl;
  
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi_v = phi;

  start_optimization = std::chrono::system_clock::now();
  Eigen::Matrix<var, Eigen::Dynamic, 1>
    theta_v = lg_logistic_solver(theta_0, phi_v, x, squared_kernel_functor(),
                                 n_samples, y);

  // compute the marginal density
  squared_kernel_functor squared_kernel;
  diff_logistic_log diff_likelihood(to_vector(n_samples), to_vector(y));

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
    Sigma_v = squared_kernel(phi_v, x, 0);
  var target = diff_likelihood.log_likelihood(theta_v);
  Eigen::Matrix<var, Eigen::Dynamic, 1> gradient, hessian;
  diff_likelihood.diff(theta_v, gradient, hessian);

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
    first_matrix = Eigen::MatrixXd::Identity(dim_theta, dim_theta)
      - diag_post_multiply(Sigma_v, hessian);

  target += - 0.5 * (log_determinant(first_matrix) 
                       + dot_product(theta_v, mdivide_left(Sigma_v, theta)));

  end_optimization = std::chrono::system_clock::now();
  elapsed_time_optimization = end_optimization - start_optimization;

  VEC g;
  AVEC parm_vec = createAVEC(phi_v(0), phi_v(1));
  target.grad(parm_vec, g);

  double obj_kinsol = diff_likelihood.log_likelihood(theta)
    - 0.5 * stan::math::quad_form_sym(Sigma, theta);
  std::cout << "LG LOGISTIC SOLVER" << std::endl
            << "obj kinsol: " << obj_kinsol << std::endl << std::endl
            << "theta: " << theta(0) << " " << theta(1) << std::endl
            << "density: " << target << std::endl
            << "autodiff grad: " << g[0] << " " << g[1] << std::endl
            << "solving and diff time: " << elapsed_time_optimization.count()
            << std::endl << std::endl;

  /////////////////////////////////////////////////////////////////////////////
  // Now do experiment using the gp_newton solver.
  start_optimization = std::chrono::system_clock::now();
  Eigen::VectorXd
    theta_gp = gp_newton_solver(theta_0, phi, x,
                                diff_logistic_log(to_vector(n_samples),
                                                  to_vector(y)),
                                squared_kernel_functor(), 1e-3, 100);
  end_optimization = std::chrono::system_clock::now();
  elapsed_time_optimization = end_optimization - start_optimization;

  double obj_gp = diff_likelihood.log_likelihood(theta_gp)
    - 0.5 * stan::math::quad_form_sym(Sigma, theta_gp);
  evaluation = f(theta_gp, phi, dat, dummy_int, 0);
  std::cout << "LG SOLVER " << std::endl
            << "target norm: " << evaluation.norm() << std::endl
            << "obj_lg: " << obj_gp << std::endl
            << "time: " << elapsed_time_optimization.count()
            << std::endl << std::endl;
  
  /////////////////////////////////////////////////////////////////////////////
  // Now do experiment using laplace_marginal function.
  using stan::math::laplace_marginal_density;

  Eigen::VectorXd theta_laplace, W_root, a, l_grad;
  Eigen::MatrixXd L;

  start_optimization = std::chrono::system_clock::now();

  double marginal_density
    = laplace_marginal_density(theta_0, phi, x,
          diff_logistic_log(to_vector(n_samples), to_vector(y)),
          squared_kernel_functor(),
          theta_laplace, W_root, L, a, l_grad,
          1e-3, 100);

  end_optimization = std::chrono::system_clock::now();
  elapsed_time_optimization = end_optimization - start_optimization;
  
  double obj_laplace = diff_likelihood.log_likelihood(theta_laplace)
    - 0.5 * stan::math::quad_form_sym(Sigma, theta_laplace);
  evaluation = f(theta_laplace, phi, dat, dummy_int, 0);
  std::cout << "LAPLACE MARGINAL SOLVER" << std::endl
            << "target norm: " << evaluation.norm() << std::endl
            << "obj_laplace: " << obj_laplace << std::endl
            << "density: " << marginal_density << std::endl
            << "time: " << elapsed_time_optimization.count()
            << std::endl << std::endl;

  /////////////////////////////////////////////////////////////////////////////
  // Throw phi as a var type
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi_v2 = phi;

  start_optimization = std::chrono::system_clock::now();
  var marginal_density_v
    = laplace_marginal_density(theta_0, phi_v2, x,
    // = laplace_marginal_density(theta, phi_v2, x,
        diff_logistic_log(to_vector(n_samples), to_vector(y)),
        squared_kernel_functor(), 1e-3, 100);

  end_optimization = std::chrono::system_clock::now();
  elapsed_time_optimization = end_optimization - start_optimization;

  VEC g2;
  AVEC parm_vec2 = createAVEC(phi_v2(0), phi_v2(1));
  marginal_density_v.grad(parm_vec2, g2);

  std::cout << "LAPLACE MARGINAL AND VARI CLASS" << std::endl
            << "density: " << value_of(marginal_density_v) << std::endl
            << "autodiff grad: " << g2[0] << " " << g2[1]
            << std::endl
            << "total time: " << elapsed_time_optimization.count()
            << std::endl << std::endl;

  // Compute derivatives using finite diff, and compare results.
  double diff = 1e-7;
  Eigen::VectorXd theta_fd = theta;  // TEST should be theta_0
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi_1u = phi_v;
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi_1l = phi_v;
  phi_1u(0) = phi_v(0) + diff;
  phi_1l(0) = phi_v(0) - diff;

  double
  marg_density_1u = laplace_marginal_density(theta_fd, value_of(phi_1u), x,
                     diff_logistic_log(to_vector(n_samples), to_vector(y)),
                     squared_kernel_functor(), 1e-3, 100);

  double
  marg_density_1l = laplace_marginal_density(theta_fd, value_of(phi_1l), x,
                     diff_logistic_log(to_vector(n_samples), to_vector(y)),
                     squared_kernel_functor(), 1e-3, 100);

  double diff_1 = (marg_density_1u - marg_density_1l) / (2 * diff);

  Eigen::Matrix<var, Eigen::Dynamic, 1> phi_2u = phi_v;
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi_2l = phi_v;
  phi_2u(1) = phi_v(1) + diff;
  phi_2l(1) = phi_v(1) - diff;

  double
  marg_density_2u = laplace_marginal_density(theta_fd, value_of(phi_2u), x,
    diff_logistic_log(to_vector(n_samples), to_vector(y)),
    squared_kernel_functor(), 1e-3, 100);

  double
  marg_density_2l = laplace_marginal_density(theta_fd, value_of(phi_2l), x,
    diff_logistic_log(to_vector(n_samples), to_vector(y)),
    squared_kernel_functor(), 1e-3, 100);

  double diff_2 = (marg_density_2u - marg_density_2l) / (2 * diff);
  std::cout << "finite diff: " << diff_1 << " " << diff_2 << std::endl;
}
