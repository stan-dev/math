#include <stan/math.hpp>

#include <stan/math/laplace/lgp_dense_system.hpp>
#include <stan/math/laplace/lgp_dense_newton_solver.hpp>
#include <stan/math/laplace/lgp_solver.hpp>
#include <stan/math/laplace/lgp_density.hpp>

#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/laplace/lgp_utility.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>


TEST(laplace, lgp_performance_density) {
  using stan::math::var;
  using stan::math::lgp_dense_newton_solver;
  using stan::math::lgp_solver;
  using stan::math::to_vector;
  using stan::math::to_row_vector;
  using stan::math::elt_multiply;
  using stan::math::multiply;
  using stan::math::transpose;
  using stan::math::exp;
  using stan::math::log_determinant;
  using stan::math::diag_post_multiply;
  using stan::math::value_of;
  using stan::math::dot_product;
  using stan::math::sum;

  std::string data_directory = "test/unit/math/laplace/data_cpp/";
  // std::string data_directory = "test/unit/math/laplace/big_data_cpp/";
  // std::string data_directory
  //   = "test/unit/math/rev/mat/functor/performance/data_cpp/";
  int n_dimensions = 5;
  Eigen::VectorXd dimensions(n_dimensions);
  dimensions << 10, 20, 50, 100, 500;

  auto start = std::chrono::system_clock::now();
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time_evaluation;
  std::chrono::duration<double> elapsed_time_jacobian;

  for (int k = 0; k < 5; k++) {
    int dim_theta = dimensions(k);

    int N = 3 * dim_theta;  // number of observations
    std::vector<int> y(N);
    std::vector<int> index(N);
    std::vector<int> sums(dim_theta);
    std::vector<int> n_samples(dim_theta);

    read_in_data(dim_theta, N, data_directory, y, index, sums, n_samples);

    int dim_phi = 2;
    Eigen::Matrix<var, Eigen::Dynamic, 1> phi(dim_phi);
    phi << 0.5, 0.9;

    start = std::chrono::system_clock::now();
    var target = 1;

    // prior on phi
    target += lognormal_lpdf(phi(0), -0.5, 0.5);
    target += normal_lpdf(phi(1), 0.5, 0.1);

    Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);
    Eigen::Matrix<var, Eigen::Dynamic, 1> theta;

    // options:
    //  1. kinsol solver
    //  2. custom Newton solver.
    //  3. algoritm 3.1 solver
    int solving_method = 1;
    auto start_optimization = std::chrono::system_clock::now();
    if (solving_method == 1)
      theta = lgp_solver(theta_0, phi, n_samples, sums, 1e-6, 100);
    if (solving_method == 2) 
      theta = lgp_dense_newton_solver(theta_0, phi, n_samples, sums, 1e-6,
                                      100, 0, 0, 1);
    auto end_optimization = std::chrono::system_clock::now();
    std::chrono::duration<double>
      elapsed_time_optimization = end_optimization - start_optimization;

    // likelihood of y conditional on phi and theta
    // (commented out: likelihood density is not relevant for 
    // comparing algorithms, since it is the same for all methods.)
    // for (int i = 0; i < N; i++)
    //   target += poisson_log_lpmf(y[i], theta(index[i] - 1));

    // ratio of conditionals on theta
    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
      Sigma = stan::math::lgp_covariance(phi, dim_theta, true);
    Eigen::Matrix<var, Eigen::Dynamic, 1> first_term(dim_theta);
    for (int i = 0; i < dim_theta; i++)  // CHECK -- use dot_product?
      first_term(i) = - n_samples[i] * exp(theta(i));

    // options:
    //  1. efficient method, see math appendix.
    //  2. former inefficient method, which should be used as a benchmark.
    //  3. method based on algorithm 3.1. of GP for ML textbook.
    int method = 1;

    if (method == 1) {
      target += - 0.5 * dot_product(theta, mdivide_left(Sigma, theta));
      {
        Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
          negative_Sigma = - Sigma;
        Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
          first_matrix = Eigen::MatrixXd::Identity(dim_theta, dim_theta) -
            diag_post_multiply(Sigma, first_term);

        target += - 0.5 * log_determinant(first_matrix);
      }
    }

    if (method == 2) {
      // Former (and much slower) method.
      // Use as a benchmark.
      Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
        Q = stan::math::inverse(Sigma);

      Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
        laplace_precision = stan::math::diag_matrix(first_term) - Q;

      target += - 0.5 * (stan::math::quad_form(Q, theta) +
        log_determinant(Sigma) +
        log_determinant(laplace_precision));
    }

    if (method == 3) {
      using stan::math::diff_poisson_log;
      var sum_log_diag_L;
      
      Eigen::Matrix<var, Eigen::Dynamic, 1>
        a = newton_update(theta, Sigma,
                          diff_poisson_log(to_vector(n_samples),
                                           to_vector(sums)),
                          sum_log_diag_L, true);

      target += -0.5 * multiply(transpose(a), theta) - sum_log_diag_L;
    }

    end = std::chrono::system_clock::now();
    elapsed_time_evaluation = end - start;

    // Compute sensitivities
    start = std::chrono::system_clock::now();
    VEC g;
    AVEC parm_vec = createAVEC(phi(0), phi(1));
    target.grad(parm_vec, g);

    end = std::chrono::system_clock::now();
    elapsed_time_jacobian = end - start;

    std::cout << "dim theta: " << dim_theta << std::endl
              << "method: " << solving_method << " " << method
              << std::endl
              << "data: " << data_directory
              << std::endl
              << "target: " << target << std::endl
              << "grad: " << g[0] << " " << g[1] << std::endl
              << "evaluation time: " << elapsed_time_evaluation.count()
              << std::endl
              << "chain time: " << elapsed_time_jacobian.count()
              << std::endl
              << "total time: " << elapsed_time_evaluation.count()
                                   + elapsed_time_jacobian.count()
              << std::endl
              << "Optimization time: " << elapsed_time_optimization.count() 
              << std::endl;

    // evaluate function with solution
    std::vector<double> dat(dim_theta * 2);
    for (int i = 0; i < dim_theta; i++) dat[i] = n_samples[i];
    for (int i = 0; i < dim_theta; i++) dat[dim_theta + i] = sums[i];
    std::vector<int> dat_int;

    inla_functor system;
    std::cout << "function eval: "
              << system(value_of(theta), value_of(phi),
                        dat, dat_int, 0).norm()
              << std::endl << std::endl;
  }
}
