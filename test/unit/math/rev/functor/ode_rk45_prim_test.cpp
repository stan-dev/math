#include <stan/math/rev.hpp>
#include <boost/numeric/odeint.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

template <typename T, stan::require_stan_scalar_t<T>* = nullptr>
T sum_(T arg) {
  return arg;
}

template <typename EigMat, stan::require_eigen_t<EigMat>* = nullptr>
auto sum_(EigMat&& arg) {
  return stan::math::sum(arg);
}

template <typename Vec, stan::require_std_vector_t<Vec>* = nullptr>
auto sum_(Vec&& arg) {
  stan::scalar_type_t<Vec> sum = 0;
  for (size_t i = 0; i < arg.size(); ++i) {
    sum += sum_(arg[i]);
  }
  return sum;
}

struct CosArg1 {
  template <typename T0, typename T1, typename... T_Args>
  inline std::vector<stan::return_type_t<T1, T_Args...>> operator()(
      const T0& t, const std::vector<T1>& y, std::ostream* msgs,
      const T_Args&... a) const {
    std::vector<typename stan::return_type<T0, T_Args...>::type> vec
        = {sum_(a)...};
    return {stan::math::cos(sum_(vec) * t)};
  }
};

struct Cos2Arg {
  template <typename T0, typename T1, typename T2, typename T3>
  inline std::vector<typename stan::return_type<T1, T2, T3>::type> operator()(
      const T0& t, const std::vector<T1>& y, std::ostream* msgs, const T2& a,
      const T3& b) const {
    return {stan::math::cos((sum_(a) + sum_(b)) * t)};
  }
};

TEST(StanMathOde_ode_rk45_tol, t0) {
  using stan::math::var;

  std::vector<double> y0 = {0.0};
  var t0 = 0.0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;

  std::vector<std::vector<var>> output = stan::math::ode_rk45_tol(
      CosArg1(), y0, t0, ts, 1e-10, 1e-10, 1e6, nullptr, a);

  output[0][0].grad();

  EXPECT_FLOAT_EQ(output[0][0].val(), 0.4165982112);
  EXPECT_FLOAT_EQ(t0.adj(), -1.0);

  stan::math::set_zero_all_adjoints();

  output[1][0].grad();

  EXPECT_FLOAT_EQ(output[1][0].val(), 0.66457668563);
  EXPECT_FLOAT_EQ(t0.adj(), -1.0);
}

TEST(StanMathOde_ode_rk45_tol, ts) {
  using stan::math::var;

  std::vector<double> y0 = {0.0};
  double t0 = 0.0;
  std::vector<var> ts = {0.45, 1.1};

  double a = 1.5;

  std::vector<std::vector<var>> output = stan::math::ode_rk45_tol(
      CosArg1(), y0, t0, ts, 1e-10, 1e-10, 1e6, nullptr, a);

  output[0][0].grad();

  EXPECT_FLOAT_EQ(output[0][0].val(), 0.4165982112);
  EXPECT_FLOAT_EQ(ts[0].adj(), 0.78070695113);

  stan::math::set_zero_all_adjoints();

  output[1][0].grad();

  EXPECT_FLOAT_EQ(output[1][0].val(), 0.66457668563);
  EXPECT_FLOAT_EQ(ts[1].adj(), -0.0791208888);
}

TEST(StanMathOde_ode_rk45_tol, scalar_arg) {
  using stan::math::var;

  std::vector<double> y0 = {0.0};
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  var a = 1.5;

  var output = stan::math::ode_rk45_tol(CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                                        nullptr, a)[0][0];

  output.grad();

  EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
  EXPECT_FLOAT_EQ(a.adj(), -0.50107310888);
}

TEST(StanMathOde_ode_rk45_tol, std_vector_arg) {
  using stan::math::var;

  std::vector<double> y0 = {0.0};
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  std::vector<var> a = {1.5};

  var output = stan::math::ode_rk45_tol(CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                                        nullptr, a)[0][0];

  output.grad();

  EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
  EXPECT_FLOAT_EQ(a[0].adj(), -0.50107310888);
}

TEST(StanMathOde_ode_rk45_tol, vector_arg) {
  using stan::math::var;

  std::vector<double> y0 = {0.0};
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  Eigen::Matrix<var, Eigen::Dynamic, 1> a(1);
  a << 1.5;

  var output = stan::math::ode_rk45_tol(CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                                        nullptr, a)[0][0];

  output.grad();

  EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
  EXPECT_FLOAT_EQ(a(0).adj(), -0.50107310888);
}

TEST(StanMathOde_ode_rk45_tol, row_vector_arg) {
  using stan::math::var;

  std::vector<double> y0 = {0.0};
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  Eigen::Matrix<var, 1, Eigen::Dynamic> a(1);
  a << 1.5;

  var output = stan::math::ode_rk45_tol(CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                                        nullptr, a)[0][0];

  output.grad();

  EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
  EXPECT_FLOAT_EQ(a(0).adj(), -0.50107310888);
}

TEST(StanMathOde_ode_rk45_tol, matrix_arg) {
  using stan::math::var;

  std::vector<double> y0 = {0.0};
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> a(1, 1);
  a << 1.5;

  var output = stan::math::ode_rk45_tol(CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                                        nullptr, a)[0][0];

  output.grad();

  EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
  EXPECT_FLOAT_EQ(a(0, 0).adj(), -0.50107310888);
}

TEST(StanMathOde_ode_rk45_tol, scalar_std_vector_args) {
  using stan::math::var;

  std::vector<double> y0 = {0.0};
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  var a0 = 0.0;
  std::vector<var> a1 = {1.5};

  var output = stan::math::ode_rk45_tol(Cos2Arg(), y0, t0, ts, 1e-8, 1e-10, 1e6,
                                        nullptr, a0, a1)[0][0];

  output.grad();

  EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
  EXPECT_FLOAT_EQ(a1[0].adj(), -0.50107310888);
}
