
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/algebra_system.hpp>
#include <stan/math/rev/mat/functor/algebra_solver_newton.hpp>
#include <stan/math/rev/mat/functor/algebra_solver_newton_custom.hpp>
#include <stan/math/rev/mat/functor/algebra_solver.hpp>
#include <stan/math/rev/mat/functor/kinsol_data.hpp>
#include <stan/math/rev/mat/functor/kinsol_solve.hpp>
#include <stan/math/laplace/lgp_solver.hpp>
#include <stan/math/laplace/lgp_dense_newton_solver.hpp>

#include <kinsol/kinsol.h>             /* access to KINSOL func., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype */
#include <sundials/sundials_math.h>    /* access to SUNRexp               */

#include <test/unit/math/rev/mat/functor/util_algebra_solver.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>

// Accessor macro
#define Ith(v, i) NV_Ith_S(v, i - 1)

// CHECK - write code without predifining constants
#define DIM 2

// CHECK - write code without defining this structure
typedef struct {
  // int dim_theta = 2;  // CHECK - can't define integer here?
  realtype n_samples[DIM];
  realtype sums[DIM];
  realtype phi;
} *UserData;

static int algebraic_system (N_Vector theta, N_Vector f, void *user_data) {
  /* Define the system of algebraic equations to be solved.
  * u: the unknown we solve for.
  * f: the returned vector.
  * data: additional coefficients.
  */
  UserData data = (UserData)user_data;
  realtype *n_samples = data->n_samples;
  realtype *sums = data->sums;
  realtype phi = data->phi;

  realtype *theta_data = N_VGetArrayPointer_Serial(theta);
  realtype theta1 = theta_data[0];
  realtype theta2 = theta_data[1];

  realtype *f_data = N_VGetArrayPointer_Serial(f);

  f_data[0] = sums[0] - n_samples[0] * exp(theta_data[0]) - theta_data[0] / phi;
  f_data[1] = sums[1] - n_samples[1] * exp(theta_data[1]) - theta_data[1] / phi;

  return (0);
}

static void PrintOutput(N_Vector u) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf(" %8.6Lg  %8.6Lg\n", Ith(u,1), Ith(u,2));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf(" %8.6g  %8.6g\n", Ith(u,1), Ith(u,2));
#else
  printf(" %8.6g  %8.6g\n", Ith(u,1), Ith(u,2));
#endif
}

TEST(matrix, kinsol) {
  // User data
  UserData data = (UserData) malloc(sizeof *data);  // CHECK - do I need this?
  data->n_samples[0] = 5;
  data->n_samples[1] = 5;
  data->sums[0] = 3;
  data->sums[1] = 10;
  data->phi = 1;

  // create initial guess and vector to store solution
  N_Vector theta = N_VNew_Serial(2);
  realtype *theta_data = N_VGetArrayPointer_Serial(theta);
  theta_data[0] = 0.0;
  theta_data[1] = 0.0;
  int global_line_search = KIN_NONE;
  int eval_jacobian = 1;  // number of steps after which the
                          // Jacobian is re-evaluated.

  // tuning parameters for the solver
  N_Vector scaling = N_VNew_Serial(DIM);
  N_VConst_Serial(1.0, scaling);  // no scaling
  realtype f_norm_tol = 1e-5;
  realtype scaling_step_tol = 1e-5;

  void *kinsol_memory;
  kinsol_memory = KINCreate();

  int flag;
  // Set tuning parameters.
  flag = KINSetFuncNormTol(kinsol_memory, f_norm_tol);
  flag = KINSetScaledStepTol(kinsol_memory, scaling_step_tol);
  flag = KINSetMaxSetupCalls(kinsol_memory, eval_jacobian);

  // Pass the system and the data structure
  flag = KINSetUserData(kinsol_memory, data);
  flag = KINInit(kinsol_memory, algebraic_system, theta);

  // Construct linear solver.
  SUNMatrix J = SUNDenseMatrix(DIM, DIM);
  SUNLinearSolver linear_solver = SUNLinSol_Dense(theta, J);
  flag = KINSetLinearSolver(kinsol_memory, linear_solver, J);

  // Call to solver  // CHECK - why do I pass scaling twice?
  flag = KINSol (kinsol_memory, theta, global_line_search, 
                 scaling, scaling);

  printf("Solutions:\n [x1, x2] = ");
  PrintOutput(theta);

  // free memory
  N_VDestroy_Serial(theta);
  N_VDestroy_Serial(scaling);
  KINFree(&kinsol_memory);
  SUNLinSolFree(linear_solver);
  SUNMatDestroy(J);
  free(data);
}

///////////////////////////////////////////////////////////////////////////////
// Now do the tests, using a regular Stan functor as a starting point.

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
    Eigen::Matrix<scalar, Eigen::Dynamic, 1> fgrad;
    int dim_theta = 2;

    Eigen::VectorXd n_samples(dim_theta);
    n_samples(0) = dat[0];
    n_samples(1) = dat[1];

    Eigen::VectorXd sums(dim_theta);
    sums(0) = dat[2];
    sums(1) = dat[3];

    return sums - stan::math::elt_multiply(n_samples,
                                           stan::math::exp(theta))
      - theta / phi(0);
  }
};

TEST(matrix, kinsol2) {
  // Apply KINSOL solver to a functor defined a la Stan (i.e. lgp_functor).
  using stan::math::kinsol_system_data;
  using stan::math::to_array_1d;
  using stan::math::kinsol_J_f;

  int dim_theta = 2;
  Eigen::VectorXd theta_0(dim_theta);
  theta_0 << 0, 0;
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

  Eigen::VectorXd phi(1);
  phi << 1;

  // tuning parameters for the solver
  int global_line_search = KIN_NONE;
  int steps_eval_jacobian = 1;  // number of steps after which the
                                // Jacobian is re-evaluated.
  N_Vector scaling = N_VNew_Serial(dim_theta);
  N_VConst_Serial(1.0, scaling);  // no scaling
  realtype f_norm_tol = 1e-5;
  realtype scaling_step_tol = 1e-5;

  /////////////////////////////////////////////////////////////////////////////
  // Build and use KINSOL solver.
  typedef kinsol_system_data<lgp_functor, kinsol_J_f> system_data;
  system_data kinsol_data(lgp_functor(), kinsol_J_f(),
                          theta_0, phi, dat, dummy_int, 0);

  void* kinsol_memory = KINCreate(); 

  int flag;  // FIX ME -- replace this with a checkflag procedure
  flag = KINInit(kinsol_memory, &system_data::kinsol_f_system, 
                 kinsol_data.nv_x_);

  // set tuning parameters -- could construct a set-option function
  flag = KINSetFuncNormTol(kinsol_memory, f_norm_tol);
  flag = KINSetScaledStepTol(kinsol_memory, scaling_step_tol);
  flag = KINSetMaxSetupCalls(kinsol_memory, steps_eval_jacobian);

  // CHECK -- why the reinterpret_cast? How does this work?
  flag = KINSetUserData(kinsol_memory, reinterpret_cast<void*>(&kinsol_data));

  // construct Linear solver
  flag = KINSetLinearSolver(kinsol_memory, kinsol_data.LS_, kinsol_data.J_);
  flag = KINSetJacFn(kinsol_memory, &system_data::kinsol_jacobian);

  // CHECK - a better way to do this conversion.
  N_Vector theta = N_VNew_Serial(dim_theta);
  realtype* theta_data = N_VGetArrayPointer_Serial(theta);
  for (int i = 0; i < dim_theta; i++) theta_data[i] = theta_0(i);

  flag = KINSol(kinsol_memory, theta,
                global_line_search, scaling, scaling);

  // CHECK - avoid / simplifies this conversion step?
  Eigen::VectorXd theta_eigen(dim_theta);
  for (int i = 0; i < dim_theta; i++) theta_eigen(i) = theta_data[i];
  EXPECT_FLOAT_EQ(-0.388925, theta_eigen(0));
  EXPECT_FLOAT_EQ( 0.628261, theta_eigen(1));
}

TEST(matrix, kinsol3) {
  // use the kinsolve_solve function.
  using stan::math::kinsol_solve;
  using stan::math::kinsol_J_f;

  int dim_theta = 2;
  Eigen::VectorXd theta_0(dim_theta);
  theta_0 << 0, 0;
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

  Eigen::VectorXd phi(1);
  phi << 1;

  Eigen::VectorXd theta = kinsol_solve(lgp_functor(), kinsol_J_f(),
                                       theta_0, phi, dat, dummy_int);

  EXPECT_FLOAT_EQ(-0.388925, theta(0));
  EXPECT_FLOAT_EQ( 0.628261, theta(1));
}

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

TEST(matrix, kinsol4) {
  using stan::math::kinsol_solve;
  using stan::math::algebra_solver;
  using stan::math::kinsol_J_f;

  int dim_theta = 2;
  Eigen::VectorXd phi(2);
  phi << 0.5, 0.9;
  Eigen::VectorXd n_samples(dim_theta);
  n_samples << 5, 5;
  Eigen::VectorXd sums(dim_theta);
  sums << 3, 10;
  std::vector<double> dat(dim_theta * 2);
  for (int i = 0; i < dim_theta; i++) dat[i] = n_samples(i);
  for (int i = 0; i < dim_theta; i++) dat[dim_theta + i] = sums(i);
  std::vector<int> dat_int;

  Eigen::VectorXd theta_0(dim_theta);
  theta_0 << 0, 0;

  // empirically determined so that the two solvers have the same "precision".
  double tol = 1e-7;
  long int max_steps = 1e+3;

  Eigen::VectorXd 
    theta = kinsol_solve(inla_functor(), kinsol_J_f(), 
                         theta_0, phi, dat, dat_int, 0,
                         tol, max_steps);

  Eigen::VectorXd
    theta_powell = algebra_solver(inla_functor(), theta_0, phi, dat, dat_int, 0,
                                  tol, tol, max_steps);

  EXPECT_FLOAT_EQ(theta_powell(0), theta(0));
  EXPECT_FLOAT_EQ(theta_powell(1), theta(1));
  
  inla_functor system;
  std::cout << "Eval newton: " << system(theta, phi, dat, dat_int, 0).transpose() << std::endl;
  std::cout << "Eval powell: " << system(theta_powell, phi, dat, dat_int, 0).transpose() << std::endl;
  
}


TEST(matrix, kinsol5) {
  // Select wich solver you want to test. The options are:
  //  1. Powell dogleg solver
  //  2. Kinsol newton solver
  //  3. Custom newton solver
  //  4. lgp solver using kinsol
  //  5. lgp solver using custom
  std::vector<bool> evaluate_solver = {0, 0, 0, 1, 1};

  // using stan::math::algebra_solver_newton;
  using stan::math::kinsol_solve;
  using stan::math::kinsol_J_f;
  using stan::math::algebra_solver_newton_custom;
  using stan::math::algebra_solver;

  // variables to monitor runtime
  auto start = std::chrono::system_clock::now();
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds_total;

  // tuning parameters
  // empirically determined so that the two solvers have the same "precision".
  double rel_tol = 1e-10;
  double fun_tol = 1e-8;
  long int max_steps = 1e+3;

  int dim_theta = 500;  // options: 10, 20, 50, 100, 500
  std::cout << "dim theta: " << dim_theta << std::endl;
  Eigen::VectorXd phi(2);
  phi << 0.5, 0.9;

  std::string data_directory = "test/unit/math/rev/mat/functor/performance/";
  std::vector<double> n_samples(dim_theta);
  std::vector<double> sums(dim_theta);

  // read data from csv file to work out n_samples and sums.
  std::ifstream input_data;
  std::string dim_theta_string = std::to_string(dim_theta);
  std::string file_m = data_directory + "data_cpp/m_" +
    dim_theta_string + ".csv";
  std::string file_sums = data_directory + "data_cpp/sums_" +
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

  std::vector<double> dat(dim_theta * 2);
  for (int i = 0; i < dim_theta; i++) dat[i] = n_samples[i];
  for (int i = 0; i < dim_theta; i++) dat[dim_theta + i] = sums[i];
  std::vector<int> dat_int;

  Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);

  // Powell solver
  Eigen::VectorXd theta;
  if (evaluate_solver[0]) {
    start = std::chrono::system_clock::now();
    theta = algebra_solver(inla_functor(), theta_0, phi, dat, dat_int);
    end = std::chrono::system_clock::now();
    elapsed_seconds_total = end - start;
    std::cout << "Time powell: " << elapsed_seconds_total.count() << std::endl;
  }

  // Kinsol solver
  Eigen::VectorXd theta_newton;
  if (evaluate_solver[1]) {
    start = std::chrono::system_clock::now();
    theta_newton 
      = kinsol_solve(inla_functor(), kinsol_J_f(), theta_0, phi, dat, dat_int);
    end = std::chrono::system_clock::now();
    elapsed_seconds_total = end - start;
    std::cout << "Time newton kinsol: " << elapsed_seconds_total.count()
              << std::endl;
  }

  // Custom Newton solver
  Eigen::VectorXd theta_newton_custom;
  if (evaluate_solver[2]) {
    start = std::chrono::system_clock::now();
    theta_newton_custom
      = algebra_solver_newton_custom(inla_functor(), theta_0, phi, dat, dat_int);
    end = std::chrono::system_clock::now();
    elapsed_seconds_total = end - start;
    std::cout << "Time newton custom: " << elapsed_seconds_total.count()
              << std::endl;
  }

  using stan::math::lgp_solver;  // wrapper around kinsol solver.

  std::vector<int> n_samples_int(dim_theta);
  for (int i = 0; i < dim_theta; i++) n_samples_int[i] = n_samples[i];

  std::vector<int> sums_int(dim_theta);
  for (int i = 0; i < dim_theta; i++) sums_int[i] = sums[i];
  Eigen::VectorXd theta_lgp;
  if (evaluate_solver[3]) {
    start = std::chrono::system_clock::now();

    // int n_elements = dim_theta * (2 + dim_theta);
    // std::vector<double> dat(n_elements);
    // for (int i = 0; i < dim_theta; i++) dat[i] = n_samples[i];
    // for (int i = 0; i < dim_theta; i++) dat[dim_theta + i] = sums[i];
    // 
    // {
    //   std::vector<double> Q_array = to_array_1d(
    //     inverse_spd(lgp_covariance(phi, dim_theta, true))
    //   );
    //   for (int i = 0; i < dim_theta * dim_theta; i++)
    //     dat[2 * dim_theta + i] = Q_array[i];
    // }
    // // start = std::chrono::system_clock::now();
    // theta_lgp = kinsol_solve(lgp_f(), lgp_J_f(), theta_0, phi, dat, dat_int);

    theta_lgp
      = lgp_solver(theta_0, phi, n_samples_int, sums_int);

    end = std::chrono::system_clock::now();
    elapsed_seconds_total = end - start;
    std::cout << "Time lgp solver: " << elapsed_seconds_total.count()
              << std::endl;
  }

  // lgp solver using custom method
  using stan::math::lgp_dense_newton_solver;

  Eigen::VectorXd theta_lgp_custom;
  if (evaluate_solver[4]) {
    start = std::chrono::system_clock::now();
    theta_lgp_custom
      = lgp_dense_newton_solver(theta_0, phi, n_samples_int, sums_int,
                                1e-6, 100, 0, 0, 1);
    // theta_lgp_custom
    //   = lgp_solver(theta_0, phi, n_samples_int, sums_int);

    end = std::chrono::system_clock::now();
    elapsed_seconds_total = end - start;
    std::cout << "Time lgp custom solver: " << elapsed_seconds_total.count()
              << std::endl;
  }

  /////////////////////////////////////////////////////////////////////////////
  // Check the solvers found a root and converged to the right solution.

  inla_functor system;

  if (evaluate_solver[0]) 
    std::cout << "powell eval: " << system(theta, phi, dat, dat_int, 0).norm() 
    << std::endl;

  if (evaluate_solver[1])
    std::cout << "newton eval: "
              << system(theta_newton, phi, dat, dat_int, 0).norm()
              << std::endl;

  if (evaluate_solver[2])
    std::cout << "custom newton eval: "
              << system(theta_newton_custom, phi, dat, dat_int, 0).norm()
              << std::endl;

  if (evaluate_solver[3])
    std::cout << "lgp solver eval: "
              << system(theta_lgp, phi, dat, dat_int, 0).norm()
              << std::endl;

  if (evaluate_solver[4])
    std::cout << "lgp solver custom eval: "
              << system(theta_lgp_custom, phi, dat, dat_int, 0).norm()
              << std::endl;

  // for (int i = 0; i < dim_theta; i++) {
  //   if (evaluate_solver(0)) EXPECT_FLOAT_EQ(theta(i), theta_newton(i));
  //   EXPECT_FLOAT_EQ(theta_newton(i), theta_newton_custom(i));
  // }
}
