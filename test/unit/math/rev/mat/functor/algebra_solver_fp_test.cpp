#include <stan/math/rev/core.hpp>
// #include <stan/math/rev/mat/functor/algebra_solver_powell.hpp>
// #include <stan/math/rev/mat/functor/algebra_solver_newton.hpp>
#include <stan/math/rev/mat/functor/algebra_solver_fp.hpp>
#include <stan/math/prim/mat/functor/finite_diff_gradient_auto.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/functor/util_algebra_solver.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <limits>

using stan::math::KinsolFixedPointEnv;
using stan::math::FixedPointSolver;
using stan::math::FixedPointADJac;
using stan::math::algebra_solver_fp;
using stan::math::value_of;
using stan::math::var;
using stan::math::to_var;
using stan::math::finite_diff_gradient_auto;

/*
 * Solve eq
 * 
 * $x \times \exp{x} - 1 = 0$
 *
 */
struct FP_exp_func_test : public ::testing::Test {
  /*
   * RHS functor
   */
  struct FP_exp_func {
    template <typename T0, typename T1>
    inline Eigen::Matrix<typename stan::return_type<T0, T1>::type, -1, 1>
    operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
               const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
               const std::vector<double>& x_r, const std::vector<int>& x_i,
               std::ostream* pstream__) const {
      using scalar = typename stan::return_type<T0, T1>::type;
      Eigen::Matrix<scalar, Eigen::Dynamic, 1> z(1);
      z(0) = stan::math::exp(-y(0) * x(0));
      return z;
    }
  };

  FP_exp_func f;
  Eigen::VectorXd x;
  Eigen::VectorXd y;
  std::vector<double> x_r;
  std::vector<int> x_i;
  std::ostream* msgs;
  std::vector<double> u_scale;
  std::vector<double> f_scale;

  void SetUp() { stan::math::recover_memory(); }

  FP_exp_func_test()
    : f(),
      x(stan::math::to_vector(std::vector<double>{0.5})),
      y(stan::math::to_vector(std::vector<double>{1.0})),
      x_r(),
      x_i(),
      msgs(nullptr),
      u_scale{1.0},
      f_scale{1.0}
  {}

  auto fd_functor(int i) {
    auto f_fd = [this, i]
      (const Eigen::VectorXd& y_) {
      KinsolFixedPointEnv<FP_exp_func> env(f, x, y_, x_r, x_i, msgs,
                                           u_scale, f_scale);
      FixedPointSolver<KinsolFixedPointEnv<FP_exp_func>, FixedPointADJac> fp;
      double f_tol = 1.e-12;
      int max_num_steps = 100;
      return fp.solve(x, y_, env, f_tol, max_num_steps)(0);
    };
    return f_fd;
  }
};

/*
 * Solve eq
 * 
 * x2 = x1^2
 * x1^2 + x2^2 = 1
 *
 */
struct FP_2d_func_test : public ::testing::Test {
  /*
   * RHS functor
   */
  struct FP_2d_func {
    template <typename T0, typename T1>
    inline Eigen::Matrix<typename stan::return_type<T0, T1>::type, -1, 1>
    operator()(const Eigen::Matrix<T0, Eigen::Dynamic, 1>& x,
               const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
               const std::vector<double>& x_r, const std::vector<int>& x_i,
               std::ostream* pstream__) const {
      using scalar = typename stan::return_type<T0, T1>::type;
      Eigen::Matrix<scalar, Eigen::Dynamic, 1> z(2);
      z(0) = y(0) * sqrt(x(1));
      z(1) = y(1) * sqrt(y(2) - x(0) * x(0));
      return z;
    }
  };

  FP_2d_func f;
  Eigen::VectorXd x;
  Eigen::VectorXd y;
  std::vector<double> x_r;
  std::vector<int> x_i;
  std::ostream* msgs;
  std::vector<double> u_scale;
  std::vector<double> f_scale;

  void SetUp() { stan::math::recover_memory(); }

  FP_2d_func_test()
    : f(),
      x(stan::math::to_vector(std::vector<double>{0.1, 0.1})),
      y(stan::math::to_vector(std::vector<double>{1.0, 1.0, 1.0})),
      x_r(),
      x_i(),
      msgs(nullptr),
      u_scale{1.0, 1.0},
      f_scale{1.0, 1.0}
  {}

  auto fd_functor(int i) {
    auto f_fd = [this, i]
      (const Eigen::VectorXd& y_) {
      KinsolFixedPointEnv<FP_2d_func> env(f, x, y_, x_r, x_i, msgs,
                                           u_scale, f_scale);
      FixedPointSolver<KinsolFixedPointEnv<FP_2d_func>, FixedPointADJac> fp;
      double f_tol = 1.e-12;
      int max_num_steps = 100;
      return fp.solve(x, y_, env, f_tol, max_num_steps)(i);
    };
    return f_fd;
  }
};

TEST_F(FP_exp_func_test, solve) {
  KinsolFixedPointEnv<FP_exp_func> env(f, x, y, x_r, x_i, msgs,
                                       u_scale, f_scale);
  FixedPointSolver<KinsolFixedPointEnv<FP_exp_func>, FixedPointADJac> fp;
  double f_tol = 1.e-12;
  int max_num_steps = 100;
  {
    Eigen::Matrix<double, -1, 1> res = fp.solve(x, y, env, f_tol, max_num_steps); // NOLINT
    EXPECT_FLOAT_EQ(res(0), 0.567143290409);
  }

  {
    x(0) = 0.1;
    y(0) = 0.8;
    Eigen::Matrix<double, -1, 1> res = fp.solve(x, y, env, f_tol, max_num_steps); // NOLINT
    EXPECT_FLOAT_EQ(res(0), 0.612584823501);
  }
}

TEST_F(FP_2d_func_test, solve) {
  KinsolFixedPointEnv<FP_2d_func> env(f, x, y, x_r, x_i, msgs,
                                      u_scale, f_scale);
  FixedPointSolver<KinsolFixedPointEnv<FP_2d_func>, FixedPointADJac> fp;
  double f_tol = 1.e-12;
  int max_num_steps = 100;
  Eigen::Matrix<double, -1, 1> res = fp.solve(x, y, env, f_tol, max_num_steps); // NOLINT
  EXPECT_NEAR(res(0), 0.7861518, 1e-5);
  EXPECT_NEAR(res(1), 0.6180333, 1e-5);
}

TEST_F(FP_exp_func_test, gradient) {
  KinsolFixedPointEnv<FP_exp_func> env(f, x, y, x_r, x_i, msgs,
                                       u_scale, f_scale);
  FixedPointSolver<KinsolFixedPointEnv<FP_exp_func>, FixedPointADJac> fp;
  double f_tol = 1.e-12;
  int max_num_steps = 100;
  Eigen::Matrix<var, -1, 1> yp(to_var(y));

  x(0) = 0.1;
  y(0) = 0.8;
  Eigen::Matrix<var, -1, 1> x_sol = fp.solve(x, yp, env, f_tol, max_num_steps); // NOLINT
  EXPECT_FLOAT_EQ(value_of(x_sol(0)), 0.612584823501);
  stan::math::set_zero_all_adjoints();
  x_sol(0).grad();
  auto f_fd = fd_functor(0);
  double fx;
  Eigen::VectorXd grad_fx;
  finite_diff_gradient_auto(f_fd, y, fx, grad_fx);
  EXPECT_FLOAT_EQ(grad_fx(0), yp(0).adj());
}

TEST_F(FP_2d_func_test, gradient) {
  KinsolFixedPointEnv<FP_2d_func> env(f, x, y, x_r, x_i, msgs,
                                      u_scale, f_scale);
  FixedPointSolver<KinsolFixedPointEnv<FP_2d_func>, FixedPointADJac> fp;
  double f_tol = 1.e-12;
  int max_num_steps = 100;
  Eigen::Matrix<var, -1, 1> yp(to_var(y));

  Eigen::Matrix<var, -1, 1> x_sol = fp.solve(x, yp, env, f_tol, max_num_steps); // NOLINT
  EXPECT_NEAR(value_of(x_sol(0)), 0.7861518, 1e-5);
  EXPECT_NEAR(value_of(x_sol(1)), 0.6180333, 1e-5);

  double fx;
  Eigen::VectorXd grad_fx;
  for (int i = 0; i < env.N_; ++i) {
    stan::math::set_zero_all_adjoints();
    x_sol(i).grad();
    auto f_fd = fd_functor(i);
    finite_diff_gradient_auto(f_fd, y, fx, grad_fx);
    for (int j = 0; j < env.M_; ++j) {
      EXPECT_FLOAT_EQ(grad_fx(j), yp(j).adj());
    }
  }
}

TEST_F(FP_2d_func_test, gradient_with_var_init_point) {
  KinsolFixedPointEnv<FP_2d_func> env(f, x, y, x_r, x_i, msgs,
                                      u_scale, f_scale);
  FixedPointSolver<KinsolFixedPointEnv<FP_2d_func>, FixedPointADJac> fp;
  double f_tol = 1.e-12;
  int max_num_steps = 100;
  Eigen::Matrix<var, -1, 1> yp(to_var(y));
  Eigen::Matrix<var, -1, 1> xp(to_var(x));

  Eigen::Matrix<var, -1, 1> x_sol = fp.solve(xp, yp, env, f_tol, max_num_steps);
  EXPECT_NEAR(value_of(x_sol(0)), 0.7861518, 1e-5);
  EXPECT_NEAR(value_of(x_sol(1)), 0.6180333, 1e-5);

  double fx;
  Eigen::VectorXd grad_fx;
  for (int i = 0; i < env.N_; ++i) {
    stan::math::set_zero_all_adjoints();
    x_sol(i).grad();
    auto f_fd = fd_functor(i);
    finite_diff_gradient_auto(f_fd, y, fx, grad_fx);
    for (int j = 0; j < env.M_; ++j) {
      EXPECT_FLOAT_EQ(grad_fx(j), yp(j).adj());
    }
  }
}

TEST_F(FP_2d_func_test, algebra_solver_fp) {
  double f_tol = 1.e-12;
  int max_num_steps = 100;

  Eigen::Matrix<double, -1, 1> xd =
    algebra_solver_fp(f, x, y, x_r, x_i, u_scale, f_scale, 0, f_tol, max_num_steps); // NOLINT
  EXPECT_NEAR(xd(0), 0.7861518, 1e-5);
  EXPECT_NEAR(xd(1), 0.6180333, 1e-5);


  Eigen::Matrix<var, -1, 1> yp(to_var(y));
  Eigen::Matrix<var, -1, 1> xp(to_var(x));
  Eigen::Matrix<var, -1, 1> xv =
    algebra_solver_fp(f, xp, yp, x_r, x_i, u_scale, f_scale, 0, f_tol, max_num_steps); // NOLINT
  EXPECT_NEAR(value_of(xv(0)), 0.7861518, 1e-5);
  EXPECT_NEAR(value_of(xv(1)), 0.6180333, 1e-5);

  double fx;
  Eigen::VectorXd grad_fx;
  const int N = x.size();
  const int M = y.size();
  for (int i = 0; i < N; ++i) {
    stan::math::set_zero_all_adjoints();
    xv(i).grad();
    auto f_fd = fd_functor(i);
    finite_diff_gradient_auto(f_fd, y, fx, grad_fx);
    for (int j = 0; j < M; ++j) {
      EXPECT_FLOAT_EQ(grad_fx(j), yp(j).adj());
    }
  }
}

TEST_F(FP_2d_func_test, exception_handling) {
  double f_tol = 1.e-12;
  int max_num_steps = 4;

  {
    std::stringstream err_msg;
    err_msg << "algebra_solver: max number of iterations: 4 exceeded.";
    std::string msg = err_msg.str();
    EXPECT_THROW_MSG(algebra_solver_fp(f, x, y, x_r, x_i, u_scale, f_scale, 0, f_tol, max_num_steps), // NOLINT
                     std::runtime_error, msg);
  }

  {
    std::stringstream err_msg;
    err_msg << "algebra_solver: initial guess has size 0";
    std::string msg = err_msg.str();
    x = Eigen::VectorXd();
    EXPECT_THROW_MSG(algebra_solver_fp(f, x, y, x_r, x_i, u_scale, f_scale, 0, f_tol, max_num_steps), // NOLINT
                     std::invalid_argument, msg);
    x = Eigen::VectorXd(2);
    x << 0.1, 0.1;
  }

  {
    std::stringstream err_msg;
    err_msg << "algebra_solver: continuous data[1] is inf";
    std::string msg = err_msg.str();
    x_r.push_back(std::numeric_limits<double>::infinity());
    EXPECT_THROW_MSG(algebra_solver_fp(f, x, y, x_r, x_i, u_scale, f_scale, 0, f_tol, max_num_steps), // NOLINT
                     std::domain_error, msg);
    x_r.clear();
  }

  {
    std::stringstream err_msg;
    err_msg << "algebra_solver: u_scale[1] is -1, but must be >= 0";
    std::string msg = err_msg.str();
    u_scale[0] = -1.0;
    EXPECT_THROW_MSG(algebra_solver_fp(f, x, y, x_r, x_i, u_scale, f_scale, 0, f_tol, max_num_steps), // NOLINT
                     std::domain_error, msg);
    u_scale[0] = 1.0;
  }

  {
    std::stringstream err_msg;
    err_msg << "algebra_solver: f_scale[1] is -1, but must be >= 0";
    std::string msg = err_msg.str();
    f_scale[0] = -1.0;
    EXPECT_THROW_MSG(algebra_solver_fp(f, x, y, x_r, x_i, u_scale, f_scale, 0, f_tol, max_num_steps), // NOLINT
                     std::domain_error, msg);
    f_scale[0] = 1.0;
  }

  {
    std::stringstream err_msg;
    err_msg << "algebra_solver: function_tolerance is -0.1, but must be >= 0";
    std::string msg = err_msg.str();
    f_tol = -0.1;
    EXPECT_THROW_MSG(algebra_solver_fp(f, x, y, x_r, x_i, u_scale, f_scale, 0, f_tol, max_num_steps), // NOLINT
                     std::domain_error, msg);
    f_tol = 1.e-8;
  }

  {
    std::stringstream err_msg;
    err_msg << "algebra_solver: size of the algebraic system's output";
    std::string msg = err_msg.str();
    x = Eigen::VectorXd::Zero(4);
    EXPECT_THROW_MSG(algebra_solver_fp(f, x, y, x_r, x_i, u_scale, f_scale, 0, f_tol, max_num_steps), // NOLINT
                     std::invalid_argument, msg);
    x = Eigen::VectorXd(2);
    x << 0.1, 0.1;
  }
}
