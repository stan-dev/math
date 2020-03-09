#include <stan/math/rev.hpp>
#include <boost/numeric/odeint.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

template <typename T>
T sum_(T arg) {
  return arg;
}

template <typename T, int RowType, int ColType>
auto sum_(const Eigen::Matrix<T, RowType, ColType>& arg) {
  return stan::math::sum(arg);
}

template <typename T>
auto sum_(const std::vector<T>& arg) {
  stan::scalar_type_t<T> sum = 0;
  for (size_t i = 0; i < arg.size(); ++i) {
    sum += sum_(arg[i]);
  }
  return sum;
}

struct CosArg1 {
  template <typename T0, typename T1, typename... T_Args>
  inline std::vector<typename stan::return_type<T1, T_Args...>::type>
  operator()(const T0& t, const std::vector<T1>& y,
	     const T_Args&... a, std::ostream* msgs) const {

    return { stan::math::cos(sum_(std::get<0>(std::make_tuple(a...))) * t) };
  }
};

struct Cos2Arg {
  template <typename T0, typename T1, typename T2, typename T3>
  inline std::vector<typename stan::return_type<T1, T2, T3>::type>
  operator()(const T0& t, const std::vector<T1>& y,
	     const T2& a, const T3& b, std::ostream* msgs) const {
    
    return { stan::math::cos((sum_(a) + sum_(b)) * t) };
  }
};

TEST(StanMathOde_ode_bdf, scalar_arg) {
  using stan::math::var;
  
  std::vector<double> y0 = { 0.0 };
  double t0 = 0.0;
  std::vector<double> ts = { 1.1 };

  var a = 1.5;
  
  var output = stan::math::ode_bdf(CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6, nullptr, a)[0][0];

  output.grad();

  EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
  EXPECT_FLOAT_EQ(a.adj(), -0.50107310888);
}

TEST(StanMathOde_ode_bdf, std_vector_arg) {
  using stan::math::var;
  
  std::vector<double> y0 = { 0.0 };
  double t0 = 0.0;
  std::vector<double> ts = { 1.1 };

  std::vector<var> a = { 1.5 };
  
  var output = stan::math::ode_bdf(CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6, nullptr, a)[0][0];

  output.grad();

  EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
  EXPECT_FLOAT_EQ(a[0].adj(), -0.50107310888);
}

TEST(StanMathOde_ode_bdf, vector_arg) {
  using stan::math::var;
  
  std::vector<double> y0 = { 0.0 };
  double t0 = 0.0;
  std::vector<double> ts = { 1.1 };

  Eigen::Matrix<var, Eigen::Dynamic, 1> a(1);
  a << 1.5;
  
  var output = stan::math::ode_bdf(CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6, nullptr, a)[0][0];

  output.grad();

  EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
  EXPECT_FLOAT_EQ(a(0).adj(), -0.50107310888);
}

TEST(StanMathOde_ode_bdf, row_vector_arg) {
  using stan::math::var;
  
  std::vector<double> y0 = { 0.0 };
  double t0 = 0.0;
  std::vector<double> ts = { 1.1 };

  Eigen::Matrix<var, 1, Eigen::Dynamic> a(1);
  a << 1.5;
  
  var output = stan::math::ode_bdf(CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6, nullptr, a)[0][0];

  output.grad();

  EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
  EXPECT_FLOAT_EQ(a(0).adj(), -0.50107310888);
}

TEST(StanMathOde_ode_bdf, matrix_arg) {
  using stan::math::var;
  
  std::vector<double> y0 = { 0.0 };
  double t0 = 0.0;
  std::vector<double> ts = { 1.1 };

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> a(1, 1);
  a << 1.5;
  
  var output = stan::math::ode_bdf(CosArg1(), y0, t0, ts, 1e-8, 1e-10, 1e6, nullptr, a)[0][0];

  output.grad();

  EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
  EXPECT_FLOAT_EQ(a(0, 0).adj(), -0.50107310888);
}

TEST(StanMathOde_ode_bdf, scalar_std_vector_args) {
  using stan::math::var;
  
  std::vector<double> y0 = { 0.0 };
  double t0 = 0.0;
  std::vector<double> ts = { 1.1 };

  var a0 = 0.0;
  std::vector<var> a1 = { 1.5 };
  
  var output = stan::math::ode_bdf(Cos2Arg(), y0, t0, ts, 1e-8, 1e-10, 1e6, nullptr,
				   a0, a1)[0][0];

  output.grad();

  EXPECT_FLOAT_EQ(output.val(), 0.66457668563);
  EXPECT_FLOAT_EQ(a1[0].adj(), -0.50107310888);
}

