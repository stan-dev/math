#include <stan/math/rev/core.hpp>
#include <stan/math/rev/mat/functor/algebra_solver_powell.hpp>
#include <stan/math/rev/mat/functor/algebra_solver_newton.hpp>
#include <stan/math/rev/mat/functor/kinsol_fp.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/functor/util_algebra_solver.hpp>
#include <test/unit/util.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using stan::math::KinsolFixedPointEnv;
using stan::math::FixedPointSolver;

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
    std::cout << "taki test: " << y(0) << " " << x(0) << " " << z(0) << "\n";
    return z;
  }
};

struct FP_exp_func_test : public ::testing::Test {
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
      x(stan::math::to_vector(std::vector<double>{0.4})),
      y(stan::math::to_vector(std::vector<double>{1.0})),
      x_r(),
      x_i(),
      msgs(nullptr),
      u_scale{1.0},
      f_scale{1.0}
  {}
};

TEST_F(FP_exp_func_test, solve) {
  KinsolFixedPointEnv<FP_exp_func> env(f, x, y, x_r, x_i, msgs,
                                       u_scale, f_scale);
  FixedPointSolver fp;
  double f_tol = 1.e-12;
  int max_num_steps = 100;
  Eigen::Matrix<double, -1, 1> res = fp.solve(x, y, env, f_tol, max_num_steps);
}
