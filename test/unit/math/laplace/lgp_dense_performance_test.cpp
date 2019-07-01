#include <stan/math/rev/core.hpp>
#include <stan/math/laplace/lgp_dense_system.hpp>
#include <stan/math/laplace/lgp_dense_newton_solver.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

TEST(laplace, lgp_performance_dim2) {
  // Reproduce simple test where the dimension of theta is 2.
  // See the laplace/lgp_dense_system_test.cpp, see test
  // lgp_newton_solver.
  using stan::math::var;
  using stan::math::value_of;
  using stan::math::lgp_dense_system;
  using stan::math::lgp_dense_newton_solver;

  // global parameter
  int dim_phi = 2;
  Eigen::VectorXd phi(dim_phi);
  phi << 1, 0.5;  // global parameters: sigma and rho

  // elements of the algebraic equation
  int dim_theta = 2;
  Eigen::VectorXd n_samples(dim_theta);
  n_samples << 5, 5;
  Eigen::VectorXd sums(dim_theta);
  sums << 3, 10;
  lgp_dense_system<double> system(phi, n_samples, sums);
  
  // tuning parameters for the algebraic solver
  Eigen::VectorXd theta_0(dim_theta);
  double tol = 1e-6;
  int max_num_steps = 1e+6;
  bool line_search = false;
  bool print_iteration = true;

  Eigen::MatrixXd theta
    = lgp_dense_newton_solver(theta_0, system, tol, max_num_steps,
                              line_search, print_iteration);
  EXPECT_FLOAT_EQ(-0.28233352, theta(0));
  EXPECT_FLOAT_EQ(0.590499, theta(1));

  // improved initial guess
  theta_0 << -0.2, 0.5;
  theta
    = lgp_dense_newton_solver(theta_0, system, tol, max_num_steps,
                              line_search, print_iteration);
  EXPECT_FLOAT_EQ(-0.28233352, theta(0));
  EXPECT_FLOAT_EQ(0.590499, theta(1));

  // poor initial guess
  theta_0 << 10, 10;
  theta
    = lgp_dense_newton_solver(theta_0, system, tol, max_num_steps,
                              line_search, print_iteration);
  EXPECT_FLOAT_EQ(-0.28233352, theta(0));
  EXPECT_FLOAT_EQ(0.590499, theta(1));
  
  // Repeat above problem, but return iteration by reference.
  // int iteration;
  // theta = lgp_dense_newton_solver(theta_0, system, iteration);
  // std::cout << "iterations (by reference): " << iteration << std::endl;

  // Initial guess influences the number of iterations
  // Neutral guess (0, 0): 9 iterations.
  // Good guess (-0.2, 0.5): 3 iterations.
  // Poor guess (10, 10): 14 iterations
  // REMARK: the tolerance also affects the number of iterations.
}

TEST(laplace, lgp_performance) {
  // Second test:
  //  use spatial data with higher dimensions.
  //  The data is generated in make_data_hd.r,
  //  with correlation 0.9 and phi = 1.0.
  //
  // NOTE: if phi = 0, solver does NOT converge. Makes sense if
  // we look at the math.

  using stan::math::var;
  using stan::math::value_of;
  using stan::math::lgp_dense_system;
  using stan::math::lgp_dense_newton_solver;

  // global parameter
  int dim_phi = 2;
  Eigen::VectorXd phi(dim_phi);
  phi << 0.5, 0.9;
  bool space_matters = true;

  int n_dimensions = 5;
  Eigen::VectorXd dimensions(n_dimensions);
  dimensions << 10, 20, 50, 100, 500;

  // Open csv file to store results.
  std::ofstream myfile;
  myfile.open("test/unit/math/laplace/deliv/results.csv");
  int n_guesses = 5;
  for (int i = 0; i < n_guesses; i++)
    myfile << i + 1 << ", ";
  myfile << n_guesses + 1 << "\n";

  // Tuning parameters for the algebraic solver
  double tol = 1e-6;
  int max_num_steps = 1e+6;

  for (int k = 0; k < n_dimensions; k++) {
    int dim_theta = dimensions(k);
    Eigen::VectorXd n_samples(dim_theta);
    Eigen::VectorXd sums(dim_theta);

    // read data from csv file (using Mike's code as a template)
    // to work out n_samples and sums.
    std::ifstream input_data;
    std::string dim_theta_string = std::to_string(dim_theta);
    std::string file_m = "test/unit/math/laplace/data_cpp/m_" +
      dim_theta_string + ".csv";
    std::string file_sums = "test/unit/math/laplace/data_cpp/sums_" +
      dim_theta_string + ".csv";

    input_data.open(file_m);
    double buffer = 0.0;
    for (int n = 0; n < dim_theta; ++n) {
      input_data >> buffer;
      n_samples(n) = buffer;
    }
    input_data.close();

    input_data.open(file_sums);
    buffer = 0.0;
    for (int n = 0; n < dim_theta; ++n) {
      input_data >> buffer;
      sums(n) = buffer;
    }
    input_data.close();

    lgp_dense_system<double> system(phi, n_samples, sums, space_matters);

    Eigen::VectorXd n_iterations(n_guesses + 1);
    int iteration = 0;

    Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);
    Eigen::VectorXd theta;
    theta = lgp_dense_newton_comp(theta_0, system, iteration, tol,
                                  max_num_steps, false);
    n_iterations(0) = iteration;

    // Guesses are based on 2.5, 25, 50, 75, and 97.5% quantiles.
    // obtained by fitting the model when dim_theta = 50.
    Eigen::VectorXd phi_grid(n_guesses);
    phi_grid << 0.22, 0.35, 0.46, 0.60, 0.95;

    for (int i = 0; i < n_guesses; i++) {
      Eigen::VectorXd grid_point(dim_phi);
      grid_point << phi_grid(i), 0.9;
      lgp_dense_system<double> system_grid(grid_point, n_samples, sums,
                                           space_matters);

      Eigen::VectorXd theta_grid
        = lgp_dense_newton_solver(theta_0, system_grid, tol, max_num_steps,
                                  false, false);

      theta = lgp_dense_newton_comp(theta_grid, system, iteration, tol, max_num_steps,
                                    false);
      n_iterations(i + 1) = iteration;
    }

    std::cout << "Dim theta: " << dim_theta 
              << " Iterations: " << n_iterations.transpose() << std::endl;

    for (int i = 0; i < n_guesses; i++)
      myfile << n_iterations(i) << ", ";
    myfile << n_iterations(n_guesses) << "\n";
  }
}
