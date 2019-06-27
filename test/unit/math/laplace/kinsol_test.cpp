
#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/algebra_system.hpp>
#include <stan/math/rev/mat/functor/kinsol_data.hpp>
#include <stan/math/rev/mat/functor/kinsol_solve.hpp>

#include <kinsol/kinsol.h>             /* access to KINSOL func., consts. */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector       */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix       */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype */
#include <sundials/sundials_math.h>    /* access to SUNRexp               */

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

// namespace stan {
// namespace math {
// 
// /**
//  * KINSOL algebraic system data holder.
//  * Based on cvodes_ode_data.
//  * (EXPERIMENTAL)
//  * 
//  * F: structure type of the functor with the system
//  * T: type of the parameter.
//  */
// template <typename F, typename T>
// class kinsol_system_data {
//   const F& f_;
//   const Eigen::VectorXd& x_;
//   const Eigen::Matrix<T, Eigen::Dynamic, 1>& y_;
//   const Eigen::VectorXd& y_dbl_;
//   const size_t N_;
//   const std::vector<double>& dat_;
//   const std::vector<int>& dat_int_;
//   std::ostream* msgs_;
// 
//   typedef kinsol_system_data<F, T> system_data;
// 
// public:
//   // system_functor object
//   // std::vector<double> x_vec_;  // CHECK - do we need this duplicate?
//   N_Vector nv_x_;
//   SUNMatrix J_;
//   SUNLinearSolver LS_;
// 
//   /**
//    * Construct KINSOL system data object.
//    * 
//    * Arguments: f, x, y, x, x_int, msgs
//    * 
//    * NOTE: might be able to specify method for constructing J?
//    * Also want to use the pre-computed precision matirx...
//    */
//   kinsol_system_data(const F& f,
//                      const Eigen::VectorXd& x,
//                      const Eigen::Matrix<T, Eigen::Dynamic, 1>& y,
//                      const std::vector<double>& dat,
//                      const std::vector<int>& dat_int,
//                      std::ostream* msgs)
//     : f_(f), x_(x), y_(y), y_dbl_(value_of(y)),
//       dat_(dat), dat_int_(dat_int), msgs_(msgs),
//       N_(x.size()), // x_vec_(to_array_1d(x)),
//       nv_x_(N_VMake_Serial(N_, &to_array_1d(x_)[0])),
//       J_(SUNDenseMatrix(N_, N_)),
//       LS_(SUNLinSol_Dense(nv_x_, J_)) { }
// 
//   ~kinsol_system_data() {
//     SUNLinSolFree(LS_);
//     SUNMatDestroy(J_);
//   }
// 
//   /**
//    * Implements the function of type kinsol_system, which is the user-defined
//    * function passed to KINSOL.
//    * 
//    * NOTE - use f instead of ydot (i.e. what gets returned by the system
//    * function).
//    */
//   static int kinsol_f_system (N_Vector x, N_Vector f, void *user_data) {
//     const system_data* explicit_system 
//       = static_cast<const system_data*>(user_data);
//     explicit_system->f_system(NV_DATA_S(x), NV_DATA_S(f));  // CHECK
// 
//     return 0;
//   }
// 
//   /**
//    * Implements the function of type CVDlsJacFn which is the user-defined
//    * callbacks for KINSOL to calculate the jacobian of the root function.
//    * The Jacobian is stored in column major format.
//    * 
//    * CHECK -- the CVODE data equivalent has a whole lot more arguments...
//    * what are they for?
//    * REMARK - tmp1 and tmp2 are pointers to memory allocated for variables
//    * of type N_Vector which can be used by KINJacFN (the function which
//    * computes the Jacobian) as temporary storage or work space.
//    * See https://computation.llnl.gov/sites/default/files/public/kin_guide-dev.pdf,
//    * page 55.
//    */
//   static int kinsol_jacobian (N_Vector x, N_Vector f,
//                               SUNMatrix J, void *user_data,
//                               N_Vector tmp1, N_Vector tmp2) {
//     const system_data* explicit_system
//       = static_cast<const system_data*>(user_data);
//     return explicit_system->jacobian_states(NV_DATA_S(x), J);
//   }
// 
// private:
//   /**
//    * Calculates the root function, using the user-supplied functor
//    * for a given value x.
//    * 
//    * NOTE: The unknown is in the double array x[], and the the output gets
//    * stored in the array f[].
//    */
//   inline void f_system(const double x[], double f[]) const {
//     const std::vector<double> x_vec(x, x + N_);
// 
//     const std::vector<double>& f_vec
//       = to_array_1d(f_(to_vector(x_vec), y_dbl_, dat_, dat_int_, msgs_));
//     std::move(f_vec.begin(), f_vec.end(), f);
// 
//     // Draft code to the above without a conversion operation
//     // const Eigen::VectorXd x_eigen
//     //   = to_vector(std::vector(x, x + N_));  // FIX ME -- do this properly.
//     // double f_eigen[N_];  // CHECK - should it be the above expression?
//     // Eigen::Map<Eigen::VectorXd>(&f_eigen[0], N_)
//     //   = f_(x_eigen, y_dbl_, dat_, dat_int_, msgs_);
//     // std::move(f_eigen.begin(), f_eigen.end(), f);
//   }
// 
//   /**
//    * Calculate the Jacobian of the system function with respect to x.
//    */
//   inline int jacobian_states(const double x[], SUNMatrix J) const {
//     const std::vector<double> x_vec(x, x + N_);
// 
//     // CHECK - do we need to nest the call to Jacobian?
//     system_functor<F, double, double, 1> 
//       system(f_, x_, y_, dat_, dat_int_, msgs_);
//     Eigen::VectorXd fx;
//     Eigen::MatrixXd Jac;
//     jacobian(system, to_vector(x_vec), fx, Jac);
// 
//     std::vector<double> jacobian_x = std::vector<double>(N_ * N_);
//     Eigen::Map<Eigen::MatrixXd>(&jacobian_x[0], N_, N_) = Jac;
// 
//     std::move(jacobian_x.begin(), jacobian_x.end(), SM_DATA_D(J));
// 
//     return 0;
//   }
// };
// 
// }  // end namespace math
// }  // end namespace stan

TEST(matrix, kinsol2) {
  // Apply KINSOL solver to a functor defined a la Stan (i.e. lgp_functor).
  using stan::math::kinsol_system_data;
  using stan::math::to_array_1d;

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
  typedef kinsol_system_data<lgp_functor> system_data;
  system_data kinsol_data(lgp_functor(), theta_0, phi, dat, dummy_int, 0);

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

  Eigen::VectorXd theta = kinsol_solve(lgp_functor(), theta_0, phi,
                                       dat, dummy_int);

  EXPECT_FLOAT_EQ(-0.388925, theta(0));
  EXPECT_FLOAT_EQ( 0.628261, theta(1));
}

