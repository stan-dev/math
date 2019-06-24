#include <stan/math.hpp>
#include <stan/math/laplace/lgp_dense_system.hpp>
#include <stan/math/laplace/lgp_dense_newton_solver.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

TEST(laplace, lgp_performance_density) {
  using stan::math::var;
  using stan::math::lgp_dense_newton_solver;
  using stan::math::to_vector;
  using stan::math::to_row_vector;
  using stan::math::elt_multiply;
  using stan::math::multiply;
  using stan::math::transpose;
  using stan::math::exp;
  using stan::math::log_determinant;
  using stan::math::diag_post_multiply;
  using stan::math::value_of;

  std::string data_directory = "test/unit/math/laplace/data_cpp/";
  int n_dimensions = 5;
  Eigen::VectorXd dimensions(n_dimensions);
  dimensions << 10, 20, 50, 100, 500;

  auto start = std::chrono::system_clock::now();
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time_evaluation;
  std::chrono::duration<double> elapsed_time_jacobian;

  for (int k = 0; k < n_dimensions; k++) {
    int dim_theta = dimensions(k);

    int N = 3 * dim_theta;  // number of observations
    std::vector<int> y(N);
    std::vector<int> index(N);
    std::vector<int> sums(dim_theta);
    std::vector<int> n_samples(dim_theta);

    // Read in data -- CHECK: is there a cleaner way of doing this?
    std::ifstream input_data;
    std::string dim_theta_string = std::to_string(dim_theta);
    std::string file_y = data_directory + "y_" + dim_theta_string + ".csv";
    std::string file_index = data_directory + "index_" +
      dim_theta_string + ".csv";
    std::string file_m = data_directory + "m_" + dim_theta_string + ".csv";
    std::string file_sums = data_directory + "sums_" +
      dim_theta_string + ".csv";

    input_data.open(file_m);
    double buffer = 0.0;
    for (int n = 0; n < dim_theta; ++n) {
      input_data >> buffer;
      n_samples[n] = buffer;
    }
    input_data.close();

    input_data.open(file_sums);
    buffer = 0.0;
    for (int n = 0; n < dim_theta; ++n) {
      input_data >> buffer;
      sums[n] = buffer;
    }
    input_data.close();

    input_data.open(file_y);
    buffer = 0.0;
    for (int n = 0; n < N; ++n) {
      input_data >> buffer;
      y[n] = buffer;
    }
    input_data.close();

    input_data.open(file_index);
    buffer = 0.0;
    for (int n = 0; n < N; ++n) {
      input_data >> buffer;
      index[n] = buffer;
    }
    input_data.close();

    int dim_phi = 2;
    Eigen::Matrix<var, Eigen::Dynamic, 1> phi(dim_phi);
    phi << 0.5, 0.9;

    start = std::chrono::system_clock::now();
    var target = 1;

    // prior on phi
    target += lognormal_lpdf(phi(0), -0.5, 0.5);
    target += normal_lpdf(phi(1), 0.5, 0.1);

    // find mode value for theta.
    bool line_search = false;
    Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);
    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
      theta = lgp_dense_newton_solver(theta_0, phi, n_samples, sums,
                                      1e-3, 100, line_search);

    // likelihood of y conditional on phi and theta
    for (int i = 0; i < N; i++)
      target += poisson_log_lpmf(y[i], theta(index[i] - 1));

    // ratio of conditionals on theta
    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
      Sigma = stan::math::lgp_covariance(phi, dim_theta, true);
    Eigen::Matrix<var, Eigen::Dynamic, 1> first_term(dim_theta);
    for (int i = 0; i < dim_theta; i++)  // CHECK -- use dot_product?
      first_term(i) = - n_samples[i] * exp(theta(i));

    bool efficient = 1;

    if (efficient) {
      target += - 0.5 * multiply(transpose(theta),
                         mdivide_left(Sigma, theta))(0);
      // target += -0.5 * log_determinant(Sigma);
      {
        Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
          negative_Sigma = - Sigma;
        Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
          first_matrix = Eigen::MatrixXd::Identity(dim_theta, dim_theta) -
            diag_post_multiply(Sigma, first_term);

        // target += -0.5 * (log_determinant(first_matrix) -
        //   log_determinant(Sigma));
        target += - 0.5 * log_determinant(first_matrix);
      }
    }
    
    if (!efficient) {
      // former (and much slower) method
      Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
        Q = stan::math::inverse(Sigma);

      Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
        laplace_precision = stan::math::diag_matrix(first_term) - Q;

      target += - 0.5 * (stan::math::quad_form(Q, theta)(0) +
        log_determinant(Sigma) +
        log_determinant(laplace_precision));
      
      // PRINT
      std::cout << value_of(log_determinant(laplace_precision)) << std::endl;
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
              << "target: " << target << std::endl
              << "grad: " << g[0] << " " << g[1] << std::endl
              << "evaluation time: " << elapsed_time_evaluation.count()
              << std::endl
              << "jacobian time: " << elapsed_time_jacobian.count()
              << std::endl << std::endl;
  }
}
