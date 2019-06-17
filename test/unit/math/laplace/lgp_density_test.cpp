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
  using stan::math::exp;

  std::string data_directory = "test/unit/math/laplace/data_cpp/";
  int n_dimensions = 5;
  Eigen::VectorXd dimensions(n_dimensions);
  dimensions << 10, 20, 50, 100, 500;

  auto start = std::chrono::system_clock::now();
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds;

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
    var target = 0;

    // phi(0) ~ logNormal(-0.5, 0.5)
    target += lognormal_lpdf(phi(0), -0.5, 0.5);

    // phi(1) ~ normal(0.5, 0.1)
    target += normal_lpdf(phi(1), 0.5, 0.1);

    // find mode value for theta.
    Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);
    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
      theta = lgp_dense_newton_solver(theta_0, phi, n_samples, sums);

    // conditional likelihood
    for (int i = 0; i < N; i++)
      target += poisson_log_lpmf(y[i], theta(index[i] - 1));

    // third component
    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
      Sigma = stan::math::lgp_covariance(phi, dim_theta, true);
    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
      Q = stan::math::inverse_spd(Sigma);

    Eigen::Matrix<var, Eigen::Dynamic, 1> first_term(dim_theta);
    for (int i = 0; i < dim_theta; i++)  // FIX ME -- use dot_product
      first_term(i) = n_samples[i] * exp(theta(i));

    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
      laplace_precision = stan::math::diag_matrix(first_term) + Q;

    // CHECK - why does quad_form not return a scalar?
    target += - 0.5 * (stan::math::quad_form(Q, theta)(0) +
      stan::math::log_determinant(Sigma) +
      stan::math::log_determinant(laplace_precision));

    // Compute sensitivities
    VEC g;
    AVEC parm_vec = createAVEC(phi(0), phi(1));
    target.grad(parm_vec, g);

    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;

    std::cout << "dim theta: " << dim_theta << std::endl
              << "target: " << target << std::endl
              << "grad: " << g[0] << " " << g[1] << std::endl
              << "elapsed time: " << elapsed_seconds.count() << std::endl
              << std::endl;
  }
}
