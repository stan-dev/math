#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <stan/math/prim/mat/fun/head.hpp>
#include <stan/math/prim/mat/fun/tail.hpp>
#include <stan/math/rev/mat/fun/mdivide_left.hpp>
#include <stan/math/prim/mat/fun/elt_multiply.hpp>
#include <stan/math/rev/mat/fun/mdivide_left.hpp>

#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/functor/util_algebra_solver.hpp>
#include <test/unit/math/rev/mat/functor/performance/util_algebra_solver_performance.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>


TEST(MathMatrix, performance_test)  {
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

  bool measure_jacobian_time_only = 1;

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
        steady_state = algebra_solver(oneCpt_functor(),
                                      init, theta, dat, n_patients);

      // compute Jacobian
      if (measure_jacobian_time_only)
        start = std::chrono::system_clock::now();

      J.resize(2 * n, 2 + 2 * n);
      for (int i = 0; i < steady_state.size(); i++) {
        if (i > 0) set_zero_all_adjoints_nested();
        grad(steady_state(i).vi_);
        for (int k = 0; k < theta.size(); k++) J(i, k) = theta(k).adj();
      }
      recover_memory_nested();

      end = std::chrono::system_clock::now();
      elapsed_seconds = end - start;
      time = elapsed_seconds.count();  // / CLOCKS_PER_SEC;  // CPU time
      file << time << ", ";
    }
    file << "\n";
  }
}

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
    using stan::math::exp;
    using stan::math::to_vector;
    using std::vector;
    using Eigen::VectorXd;
    using Eigen::MatrixXd;
    using stan::math::head;
    using stan::math::tail;

    int n_groups = theta.size();
    VectorXd n_samples = to_vector(head(dat, n_groups));
    VectorXd sums = to_vector(tail(dat, dat.size() - n_groups));
    Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>
      Sigma = covariance(parm, n_groups, 1);

    return sums - stan::math::elt_multiply(n_samples, exp(theta)) -
      stan::math::mdivide_left(Sigma, theta);
  }
};

/*
TEST(MathMatrix, performance_test2)  {
  // simple test for INLA problem with Poisson GLM.
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

  // Open csv file to store results.
  std::string working_directory = "test/unit/math/rev/mat/functor/performance/";
  std::ofstream myfile;
  myfile.open(working_directory + "results.csv");
  int n_guesses = 5;
  for (int i = 0; i < n_guesses; i++)
    myfile << i + 1 << ", ";
  myfile << n_guesses + 1 << "\n";

  // variables to meausre time
  auto start = std::chrono::system_clock::now();
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds;

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
    Matrix<var, Dynamic, 1> 
      theta = algebra_solver(inla_functor(), theta_0, parm,
                             dat, dat_int);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    std::cout << "Evaluation time: " << elapsed_seconds.count() << std::endl;

    // compute Jacobian. Use a dummy variable to create the right
    // initial cotangent vector.
    var dummy_cot = sum(theta);
    VEC g;
    AVEC parm_vec = createAVEC(parm(0), parm(1));

    start = std::chrono::system_clock::now();
    dummy_cot.grad(parm_vec, g);
    end = std::chrono::system_clock::now();
    elapsed_seconds = end - start;
    std::cout << "Jacobian time: " << elapsed_seconds.count() << std::endl;
    
  }
}  */
