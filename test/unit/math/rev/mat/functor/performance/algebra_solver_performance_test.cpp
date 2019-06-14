#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <stan/math/rev/mat/functor/algebra_solver_newton.hpp>
#include <stan/math/rev/mat/fun/mdivide_left.hpp>
#include <stan/math/rev/mat/fun/mdivide_left.hpp>
// #include <stan/math/rev/mat/fun/multiply.hpp>
#include <stan/math/rev/mat.hpp>  // CHECK -- do I need to include this?
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <stan/math/prim/mat/fun/head.hpp>
#include <stan/math/prim/mat/fun/tail.hpp>
#include <stan/math/prim/mat/fun/elt_multiply.hpp>
#include <stan/math/laplace/lgp_dense_system.hpp>
#include <stan/math/laplace/lgp_dense_newton_solver.hpp>

#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/functor/util_algebra_solver.hpp>
#include <test/unit/math/rev/mat/functor/performance/util_algebra_solver_performance.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>

/*
TEST(MathMatrix, performance_test)  {
  // FIX ME - the test should only compute gradients once, computing
  // multiple gradients is a waste of time.

  // Additional code to measure the runtime for solving an alegbraic
  // equation and computing its Jacobian.
  // Based on the pharmacometrics experiment in the autodiff review paper.

  using stan::math::var;
  using stan::math::rep_vector;
  using stan::math::algebra_solver;
  using stan::math::start_nested;
  using stan::math::set_zero_all_adjoints_nested;
  using stan::math::recover_memory_nested;
  
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using Eigen::VectorXd;

  using std::cout;
  using std::vector;

  std::srand(1954);
  std::default_random_engine generator;
  std::uniform_real_distribution<double> unif(0.7, 1.3);

  bool measure_jacobian_time_only = 0;

  // variables to meausre time
  auto start = std::chrono::system_clock::now();
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds;
  double time;

  int N = 100;
  vector<double> alpha_dbl(2 * N);
  for (int i = 0; i < 2*N; i++) alpha_dbl[i] = unif(generator);

  // open file to save results
  std::ofstream file;
  file.open("test/unit/math/rev/mat/functor/performance/algebra_performance.csv");
  for (int i = 1; i <= 100; i++) file << i <<", ";
  file << "\n";

  for (int n = 1; n <= N; n++) {
    cout << n << " / " << N << "\n";

    // generate variables which require no sensitivities
    vector<int> n_patients = {n};
    vector<double> dat = {1000, 1}; // 1000 g, administered every hour.
    Eigen::VectorXd init(2*n);
    init = rep_vector(500, 2 * n);  // initial guess

    Matrix<double, Dynamic, Dynamic> J;
    

    for (int j = 0; j <100; j++) {
      if (!measure_jacobian_time_only) 
        start = std::chrono::system_clock::now();

      start_nested();;
      Matrix<var, Dynamic, 1> theta(2 + 2 * n);
      theta(0) = 0.5;  // k1 for population.
      theta(1) = 0.2;  // k2 for population.
      for (int l = 2; l < theta.size(); l++) theta(l) = alpha_dbl[l - 2];

      Matrix<var, Dynamic, 1>
        steady_state = algebra_solver_newton(oneCpt_functor(),
                                             init, theta, dat, n_patients);
      // Matrix<var, Dynamic, 1>
      //   steady_state = algebra_solver(oneCpt_functor(),
      //                                 init, theta, dat, n_patients);

      // compute Jacobian
      if (measure_jacobian_time_only)
        start = std::chrono::system_clock::now();
      
      // it suffices to take the gradient with respect to one output,
      // since the whole Jacobian gets computed.

      grad(steady_state(1).vi_);

      // J.resize(2 * n, 2 + 2 * n);
      // for (int i = 0; i < steady_state.size(); i++) {
      //   if (i > 0) set_zero_all_adjoints_nested();
      //   grad(steady_state(i).vi_);
      //   for (int k = 0; k < theta.size(); k++) J(i, k) = theta(k).adj();
      // }
      // recover_memory_nested();

      end = std::chrono::system_clock::now();
      elapsed_seconds = end - start;
      time = elapsed_seconds.count();  // / CLOCKS_PER_SEC;  // CPU time
      file << time << ", ";
    }
    file << "\n";
  }
} */

///////////////////////////////////////////////////////////////////////////////
/* Test solver for INLA optimization problem with Poisson GLM. */
template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
covariance (Eigen::Matrix<T, Eigen::Dynamic, 1> phi, int M,
            bool space_matters = false) {
  using std::pow;
  T sigma = phi[0];
  T rho = phi[1];
  double exponent;

  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Sigma(M, M);

  for (int i = 0; i < M; i++) {
    for (int j = 0; j < i; j++) {
      if (space_matters) {exponent = i - j;} else {exponent = 1;}
      Sigma(i, j) = pow(rho, exponent) * sigma;
      Sigma(j, i) = Sigma(i, j);
    }
    Sigma(i, i) = sigma;
  }

  return Sigma;
}


struct inla_functor {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename stan::return_type<T0, T1>::type,
                       Eigen::Dynamic, 1>
  operator() (const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
              const Eigen::Matrix<T1, Eigen::Dynamic, 1>& parm,
              const std::vector<double>& dat,
              const std::vector<int>& dat_int,
              std::ostream* pstream__ = 0) const {
    using stan::math::to_vector;
    using stan::math::head;
    using stan::math::tail;

    int n_groups = theta.size();
    Eigen::VectorXd n_samples = to_vector(head(dat, n_groups));
    Eigen::VectorXd sums = to_vector(tail(dat, dat.size() - n_groups));
    Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>
      Sigma = covariance(parm, n_groups, 1);

    return sums - stan::math::elt_multiply(n_samples, stan::math::exp(theta)) -
      stan::math::mdivide_left(Sigma, theta);
  }
};

struct inla_functor2 {
  template <typename T0, typename T1>
  inline Eigen::Matrix<typename stan::return_type<T0, T1>::type,
                       Eigen::Dynamic, 1>
  operator() (const Eigen::Matrix<T0, Eigen::Dynamic, 1>& theta,
           const Eigen::Matrix<T1, Eigen::Dynamic, 1>& parm,
           const std::vector<double>& dat,
           const std::vector<int>& dat_int,
           std::ostream* pstream__ = 0) const {
    using stan::math::to_vector;
    using stan::math::head;
    using stan::math::tail;

    int n_groups = theta.size();
    Eigen::VectorXd n_samples = to_vector(head(dat, n_groups));
    Eigen::VectorXd sums = to_vector(tail(dat, dat.size() - n_groups));

    Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>
      sigma_inv(n_groups, n_groups);
    int index_offset = 0;
    for (int i = 0; i < n_groups; i++) {
      index_offset += i; 
      for (int j = 0; j <= i; j++) {
        sigma_inv(i, j) = parm(index_offset + j);
        sigma_inv(j, i) = sigma_inv(i, j);
      }
    }

    return sums - stan::math::elt_multiply(n_samples, stan::math::exp(theta)) -
      stan::math::multiply(sigma_inv, theta);
  }
};
/*
TEST(MathMatrix, optimization_inla) {
  using std::vector;
  using stan::math::var;

  int dim_theta = 2;
  Eigen::VectorXd theta(dim_theta);
  theta << -0.7268203, 1.3347728;

  Eigen::VectorXd phi(2);
  phi << 4, 0;

  Eigen::VectorXd n_samples(dim_theta);
  n_samples << 5, 5;
  Eigen::VectorXd sums(dim_theta);
  sums << 3, 10;

  // Test evaluation of gradient of the log density, when using
  // inla_functor.
  // TO DO: do the same test on inla_functor2
  inla_functor system;
  vector<double> dat(dim_theta * 2);
  for (int i = 0; i < dim_theta; i++) dat[i] = n_samples(i);
  for (int i = 0; i < dim_theta; i++) dat[dim_theta + i] = sums(i);
  vector<int> dat_int;

  EXPECT_FLOAT_EQ(0.7644863, system(theta, phi, dat, dat_int)(0));
  EXPECT_FLOAT_EQ(-9.3293567, system(theta, phi, dat, dat_int)(1));

  // Compute the solution and the sensitivities using Powell's
  // method and use these results as a benchmark.
  Eigen::MatrixXd solver_gradient(2, 2);
  Eigen::VectorXd powell_solution;
  Eigen::VectorXd theta_0(dim_theta);  // initial guess
  theta_0 << 0, 0;

  for (int k = 0; k < dim_theta; k++) {
    var sigma = 0.5;
    var corr = 0.9;
    Eigen::Matrix<var, Eigen::Dynamic, 1> phi_v(2);
    phi_v << sigma, corr;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
      = algebra_solver(inla_functor(), theta_0, phi_v, dat, dat_int);
    powell_solution = value_of(theta);  // a bit redundant...

    AVEC parameters = createAVEC(phi_v(0), phi_v(1));
    VEC g;
    theta(k).grad(parameters, g);
    solver_gradient(k, 0) = g[0];
    solver_gradient(k, 1) = g[1];
  }
  
  ////////////////////////////////////////////////////////////////////////
  // TEST: use Newton's solver and the inla_functor.
  for (int k = 0; k < dim_theta; k++) {
    var sigma = 0.5;
    var corr = 0.9;
    Eigen::Matrix<var, Eigen::Dynamic, 1> phi_v(2);
    phi_v << sigma, corr;

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
      = algebra_solver_newton(inla_functor(), theta_0, phi_v, dat, dat_int);

    AVEC parameters = createAVEC(phi_v(0), phi_v(1));
    VEC g;
    theta(k).grad(parameters, g);

    EXPECT_FLOAT_EQ(powell_solution(0), theta(0).val());
    EXPECT_FLOAT_EQ(powell_solution(1), theta(1).val());
    EXPECT_FLOAT_EQ(solver_gradient(k, 0), g[0]);
    EXPECT_FLOAT_EQ(solver_gradient(k, 1), g[1]);
  }

  ////////////////////////////////////////////////////////////////////////
  // TEST: use the newton solver, this time with inla_functor2.
  // In this configuration, the precision matrix (inverse covariance
  // matrix) is pre-computed, to avoid doing an mdivide operation
  // inside the algebraic operation.
  inla_functor2 system2;
  int triangle_size = 0.5 * dim_theta * (dim_theta + 1);

  // make sure we recover the solution stored in Powell solution
  for (int k = 0; k < dim_theta; k++) {
    var sigma = 0.5;
    var corr = 0.9;
    Eigen::Matrix<var, Eigen::Dynamic, 1> phi_v(2);
    phi_v << sigma, corr;
    Eigen::Matrix<var, Eigen::Dynamic, 1> parm(triangle_size);
    Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
      sigma_inv = stan::math::inverse_spd(covariance(phi_v, dim_theta, 1));
    int index_offset = 0;
    for (int i = 0; i <= dim_theta - 1; i++) {
      index_offset += i;
      for (int j = 0; j <= i; j++)
        parm(index_offset + j) = sigma_inv(i, j);
    }

    Eigen::Matrix<var, Eigen::Dynamic, 1> theta
      = algebra_solver_newton(inla_functor2(), theta_0, parm, dat, dat_int);

    AVEC parameters = createAVEC(phi_v(0), phi_v(1));
    VEC g;
    theta(k).grad(parameters, g);
    EXPECT_FLOAT_EQ(powell_solution(0), theta(0).val());
    EXPECT_FLOAT_EQ(powell_solution(1), theta(1).val());
    EXPECT_FLOAT_EQ(solver_gradient(k, 0), g[0]);
    EXPECT_FLOAT_EQ(solver_gradient(k, 1), g[1]);
  }

  ////////////////////////////////////////////////////////////////////////
  // TEST: Use lgp_newton solver function.
  double tol = 1e-6;
  int max_num_steps = 1e+3;
  bool line_search = false;

  for (int k = 0; k < dim_theta; k++) {
    var sigma = 0.5;
    var corr = 0.9;
    Eigen::Matrix<var, Eigen::Dynamic, 1> phi_v(2);
    phi_v << sigma, corr;

    stan::math::lgp_dense_system<double> system_dbl(value_of(phi_v), n_samples, sums,
                                               true);

    Eigen::Matrix<var, Eigen::Dynamic, 1>
      theta = lgp_dense_newton_solver(theta_0, phi_v, system_dbl);

    AVEC parameters = createAVEC(phi_v(0), phi_v(1));
    VEC g;
    theta(k).grad(parameters, g);
    EXPECT_FLOAT_EQ(powell_solution(0), theta(0).val());
    EXPECT_FLOAT_EQ(powell_solution(1), theta(1).val());
    EXPECT_FLOAT_EQ(solver_gradient(k, 0), g[0]);
    EXPECT_FLOAT_EQ(solver_gradient(k, 1), g[1]);
  }
} */

/*
TEST(MathMatrix, performance_test_inla)  {
  // performance test for INLA problem with Poisson GLM.
  using stan::math::var;
  using stan::math::append_row;
  using std::vector;
  using Eigen::VectorXd;
  using Eigen::MatrixXd;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::head;
  using stan::math::tail;

  std::string solver = "powell";

  int dim_phi = 2;
  Eigen::VectorXd phi(dim_phi);
  phi << 0.5, 0.9;
  bool space_matters = true;

  int n_dimensions = 5;
  Eigen::VectorXd dimensions(n_dimensions);
  // dimensions << 10, 20, 50, 100, 500;
  dimensions << 10, 20, 20, 20, 500;
  
  // tuning parameter for the Powell solver
  double rel_tol = 1e-12;
  double ftol = 1e-6;
  int max_num_steps = 1e4;

  std::string working_directory = "test/unit/math/rev/mat/functor/performance/";

  // variables to meausre time
  auto start = std::chrono::system_clock::now();
  auto start_jacobian = std::chrono::system_clock::now();
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds_total;
  std::chrono::duration<double> elapsed_seconds_jacobian;

  // Note: the last tier, with dim = 500 takes 20 minutes to run,
  // so it's best to leave that out.
  for (int k = 0; k < n_dimensions; k++) {

    std::cout << "Dimension: " << dimensions(k) << std::endl;

    int dim_theta = dimensions(k);
    Eigen::VectorXd n_samples(dim_theta);
    Eigen::VectorXd sums(dim_theta);

    // read data from csv file to work out n_samples and sums.
    std::ifstream input_data;
    std::string dim_theta_string = std::to_string(dim_theta);
    std::string file_m = working_directory + "data_cpp/m_" +
      dim_theta_string + ".csv";
    std::string file_sums = working_directory + "data_cpp/sums_" +
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

    // solve algebraic equation
    VectorXd theta_0 = VectorXd::Zero(dim_theta);
    Matrix<var, Dynamic, 1> parm = phi;

    vector<double> dat(dim_theta * 2);
    for (int i = 0; i < dim_theta; i++) dat[i] = n_samples(i);
    for (int i = 0; i < dim_theta; i++) dat[dim_theta + i] = sums(i);
    vector<int> dat_int;
    
    start = std::chrono::system_clock::now();

    Matrix<var, Dynamic, 1> theta;
    if (solver == "newton") {
      theta = algebra_solver_newton(inla_functor(), theta_0, parm,
                                    dat, dat_int);
    } else {
      theta = algebra_solver(inla_functor(), theta_0, parm,
                             dat, dat_int, 0, rel_tol, ftol, max_num_steps);
    }

    end = std::chrono::system_clock::now();
    elapsed_seconds_total = end - start;
    elapsed_seconds_jacobian = end - start;

    // compute Jacobian. Use a dummy variable to create the right
    // initial cotangent vector.
    var dummy_cot = sum(theta);
    VEC g;
    AVEC parm_vec = createAVEC(parm(0), parm(1));

    start_jacobian = std::chrono::system_clock::now();
    dummy_cot.grad(parm_vec, g);
    end = std::chrono::system_clock::now();

    elapsed_seconds_total = end - start;
    elapsed_seconds_jacobian = end - start_jacobian; 
    std::cout << "Solver: " << solver << std::endl
              << "Total time: " << elapsed_seconds_total.count() 
              << std::endl
              << "Jacobian time: "
              << elapsed_seconds_jacobian.count()
              << std::endl;
    
    std::cout << "theta: " << std::endl;
    for (int i = 0; i < theta.size(); i++) std::cout << theta(i).val() << " ";
    std::cout << std::endl;
    for (size_t i = 0; i < g.size(); i++) std::cout << g[i] << " ";
    std::cout << std::endl << std::endl;
  }
}  */

/*
TEST(MathMatrix, performance_test_inla2)  {
  // simple test for INLA problem with Poisson GLM.
  // This time precompute the precision matrix.
  using stan::math::var;
  using stan::math::algebra_solver_newton;
  using stan::math::append_row;
  using std::vector;
  using Eigen::VectorXd;
  using Eigen::MatrixXd;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::head;
  using stan::math::tail;

  bool space_matters = true;
  int dim_phi = 2;

  int n_dimensions = 5;
  Eigen::VectorXd dimensions(n_dimensions);
  dimensions << 10, 20, 50, 100, 500;

  // Open csv file to store results.
  std::string working_directory = "test/unit/math/rev/mat/functor/performance/";

  // variables to meausre time
  auto start = std::chrono::system_clock::now();
  auto start_jacobian = std::chrono::system_clock::now();
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds_total;
  std::chrono::duration<double> elapsed_seconds_jacobian;

  // Note: the last tier, with dim = 500 takes 20 minutes to run,
  // so it's best to leave that out.
  for (int k = 0; k < n_dimensions; k++) {

    std::cout << "Dimension: " << dimensions(k) << std::endl;

    int dim_theta = dimensions(k);
    Eigen::VectorXd n_samples(dim_theta);
    Eigen::VectorXd sums(dim_theta);

    // read data from csv file to work out n_samples and sums.
    std::ifstream input_data;
    std::string dim_theta_string = std::to_string(dim_theta);
    std::string file_m = working_directory + "data_cpp/m_" +
      dim_theta_string + ".csv";
    std::string file_sums = working_directory + "data_cpp/sums_" +
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

    // solve algebraic equation
    VectorXd theta_0 = VectorXd::Zero(dim_theta);

    vector<double> dat(dim_theta * 2);
    for (int i = 0; i < dim_theta; i++) dat[i] = n_samples(i);
    for (int i = 0; i < dim_theta; i++) dat[dim_theta + i] = sums(i);
    vector<int> dat_int;

    start = std::chrono::system_clock::now();

    Eigen::Matrix<var, Eigen::Dynamic, 1> phi(dim_phi);
    phi << 0.5, 0.9;
    int triangle_size = 0.5 * dim_theta * (dim_theta + 1);
    Matrix<var, Dynamic, 1> parm(triangle_size);
    Matrix<var, Dynamic, Dynamic> 
      sigma_inv = stan::math::inverse_spd(covariance(phi, dim_theta, 1));
    int index_offset = 0;
    for (int i = 0; i <= dim_theta - 1; i++) {
      index_offset += i;
      for (int j = 0; j <= i; j++)
        parm(index_offset + j) = sigma_inv(i, j);
    }

    Matrix<var, Dynamic, 1>
      theta = algebra_solver_newton(inla_functor2(), theta_0, parm,
                                    dat, dat_int);

    end = std::chrono::system_clock::now();
    elapsed_seconds_total = end - start;

    // compute Jacobian. Use a dummy variable to create the right
    // initial cotangent vector.
    var dummy_cot = sum(theta);
    VEC g;
    AVEC parm_vec = createAVEC(phi(0), phi(1));

    start_jacobian = std::chrono::system_clock::now();
    dummy_cot.grad(parm_vec, g);
    end = std::chrono::system_clock::now();

    elapsed_seconds_total = end - start;
    elapsed_seconds_jacobian = end - start_jacobian; 
    std::cout << "Total time: " << elapsed_seconds_total.count() 
              << std::endl
              << "Jacobian time: "
              << elapsed_seconds_jacobian.count()
              << std::endl;

    std::cout << "theta: " << std::endl;
    for (int i = 0; i < theta.size(); i++) std::cout << theta(i).val() << " ";
    std::cout << std::endl;
    for (size_t i = 0; i < g.size(); i++) std::cout << g[i] << " ";
    std::cout << std::endl << std::endl;
  }
} */

TEST(MathMatrix, performance_test_inla3) {
  // Use lgp_dense_newton_solver, which offers analytical
  // derivatives (check!)
  using stan::math::var;
  using stan::math::append_row;
  using std::vector;
  using Eigen::VectorXd;
  using Eigen::MatrixXd;
  using Eigen::Matrix;
  using Eigen::Dynamic;
  using stan::math::head;
  using stan::math::tail;

  int dim_phi = 2;
  Eigen::VectorXd phi(dim_phi);
  phi << 0.5, 0.9;
  bool space_matters = true;
  
  int n_dimensions = 5;
  Eigen::VectorXd dimensions(n_dimensions);
  dimensions << 10, 20, 50, 100, 500;
  
  // tuning parameter for the solver
  double tol = 1e-6;
  int max_num_steps = 1e+3;

  std::string working_directory = "test/unit/math/rev/mat/functor/performance/";

  // variables to meausre time
  auto start = std::chrono::system_clock::now();
  auto start_jacobian = std::chrono::system_clock::now();
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds_total;
  std::chrono::duration<double> elapsed_seconds_jacobian;

  for (int k = 0; k < n_dimensions; k++) {
    std::cout << "Dimension: " << dimensions(k) << std::endl;

    int dim_theta = dimensions(k);
    Eigen::VectorXd n_samples(dim_theta);
    Eigen::VectorXd sums(dim_theta);

    // read data from csv file to work out n_samples and sums.
    std::ifstream input_data;
    std::string dim_theta_string = std::to_string(dim_theta);
    std::string file_m = working_directory + "data_cpp/m_" +
      dim_theta_string + ".csv";
    std::string file_sums = working_directory + "data_cpp/sums_" +
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

    // solve algebraic equation
    VectorXd theta_0 = VectorXd::Zero(dim_theta);
    Matrix<var, Dynamic, 1> parm = phi;

    vector<double> dat(dim_theta * 2);
    for (int i = 0; i < dim_theta; i++) dat[i] = n_samples(i);
    for (int i = 0; i < dim_theta; i++) dat[dim_theta + i] = sums(i);
    vector<int> dat_int;
    
    start = std::chrono::system_clock::now();
    stan::math::lgp_dense_system<double> system_dbl(value_of(parm), n_samples, sums,
                                               space_matters = true);

    // Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>
    //   Q_ = stan::math::inverse_spd(stan::math::lgp_covariance(phi, 500, true));

    // start = std::chrono::system_clock::now();
    Eigen::Matrix<var, Eigen::Dynamic, 1> // theta = theta_0;
      theta = lgp_dense_newton_solver(theta_0, parm, system_dbl);

    end = std::chrono::system_clock::now();
    elapsed_seconds_total = end - start;

    // compute Jacobian. Use a dummy variable to create the right
    // initial cotangent vector.
    var dummy_cot = sum(theta);
    VEC g;
    AVEC parm_vec = createAVEC(parm(0), parm(1));

    start_jacobian = std::chrono::system_clock::now();
    dummy_cot.grad(parm_vec, g);
    end = std::chrono::system_clock::now();
    
    elapsed_seconds_total = end - start;
    elapsed_seconds_jacobian = end - start_jacobian; 
    std::cout << "Total time: " << elapsed_seconds_total.count() 
              << std::endl
              << "Jacobian time: "
              << elapsed_seconds_jacobian.count()
              << std::endl;

    std::cout << "theta: " << std::endl;
    for (int i = 0; i < theta.size(); i++) std::cout << theta(i).val() << " ";
    std::cout << std::endl;
    for (size_t i = 0; i < g.size(); i++) std::cout << g[i] << " ";
    std::cout << std::endl << std::endl;
  }
}
