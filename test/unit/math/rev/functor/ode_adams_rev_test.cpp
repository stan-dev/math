#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>

using stan::math::var;

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
  template <typename T0, typename T_y, typename... T_Args>
  inline Eigen::Matrix<stan::return_type_t<T_y, T_Args...>, Eigen::Dynamic, 1>
  operator()(const T0& t, const T_y& y, std::ostream* msgs,
             const T_Args&... a) const {
    std::vector<typename stan::return_type<T0, T_Args...>::type> vec
        = {sum_(a)...};
    Eigen::Matrix<stan::return_type_t<T_y, T_Args...>, Eigen::Dynamic, 1> out(
        1);
    out << stan::math::cos(sum_(vec) * t);
    return out;
  }
};

struct Cos2Arg {
  template <typename T0, typename T_y, typename T2, typename T3>
  inline Eigen::Matrix<stan::return_type_t<T_y, T2, T3>, Eigen::Dynamic, 1>
  operator()(const T0& t, const Eigen::Matrix<T_y, Eigen::Dynamic, 1>& y,
             std::ostream* msgs, const T2& a, const T3& b) const {
    Eigen::Matrix<stan::return_type_t<T_y, T2, T3>, Eigen::Dynamic, 1> out(1);
    out << stan::math::cos((sum_(a) + sum_(b)) * t);
    return out;
  }
};

TEST(StanMathOde_ode_adams, int_t0) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  int t0 = 0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> output
      = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a);

  EXPECT_NEAR(output[0][0], 0.4165982112, 1e-5);
  EXPECT_NEAR(output[1][0], 0.66457668563, 1e-5);
}

TEST(StanMathOde_ode_adams, int_ts) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0.0;
  std::vector<int> ts = {1, 2};

  double a = 1.5;

  std::vector<Eigen::Matrix<double, Eigen::Dynamic, 1>> output
      = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a);

  EXPECT_NEAR(output[0][0], 0.6649966577, 1e-5);
  EXPECT_NEAR(output[1][0], 0.09408000537, 1e-5);
}

TEST(StanMathOde_ode_adams, t0) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  var t0 = 0.0;
  std::vector<double> ts = {0.45, 1.1};

  double a = 1.5;

  std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> output
      = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a);

  output[0][0].grad();

  EXPECT_NEAR(output[0][0].val(), 0.4165982112, 1e-5);
  EXPECT_NEAR(t0.adj(), -1.0, 1e-5);

  stan::math::set_zero_all_adjoints();

  output[1][0].grad();

  EXPECT_NEAR(output[1][0].val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(t0.adj(), -1.0, 1e-5);
}

TEST(StanMathOde_ode_adams, ts) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0.0;
  std::vector<var> ts = {0.45, 1.1};

  double a = 1.5;

  std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> output
      = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a);

  output[0][0].grad();

  EXPECT_NEAR(output[0][0].val(), 0.4165982112, 1e-5);
  EXPECT_NEAR(ts[0].adj(), 0.78070695113, 1e-5);

  stan::math::set_zero_all_adjoints();

  output[1][0].grad();

  EXPECT_NEAR(output[1][0].val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(ts[1].adj(), -0.0791208888, 1e-5);
}

TEST(StanMathOde_ode_adams, ts_repeat) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0.0;
  std::vector<var> ts = {0.45, 0.45, 1.1, 1.1};

  double a = 1.5;

  std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> output
      = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a);

  EXPECT_EQ(output.size(), ts.size());

  output[0][0].grad();

  EXPECT_NEAR(output[0][0].val(), 0.4165982112, 1e-5);
  EXPECT_NEAR(ts[0].adj(), 0.78070695113, 1e-5);

  stan::math::set_zero_all_adjoints();

  output[1][0].grad();

  EXPECT_NEAR(output[1][0].val(), 0.4165982112, 1e-5);
  EXPECT_NEAR(ts[1].adj(), 0.78070695113, 1e-5);

  stan::math::set_zero_all_adjoints();

  output[2][0].grad();

  EXPECT_NEAR(output[2][0].val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(ts[2].adj(), -0.0791208888, 1e-5);

  stan::math::set_zero_all_adjoints();

  output[3][0].grad();

  EXPECT_NEAR(output[3][0].val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(ts[3].adj(), -0.0791208888, 1e-5);
}

TEST(StanMathOde_ode_adams, scalar_arg) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  var a = 1.5;

  var output = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a)[0][0];

  output.grad();

  EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(a.adj(), -0.50107310888, 1e-5);
}

TEST(StanMathOde_ode_adams, scalar_arg_multi_time) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0.0;
  std::vector<double> ts = {0.45, 1.1};

  var a = 1.5;

  std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> output
      = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a);

  output[0](0).grad();

  EXPECT_NEAR(output[0](0).val(), 0.4165982112, 1e-5);
  EXPECT_NEAR(a.adj(), -0.04352005542, 1e-5);

  stan::math::set_zero_all_adjoints();

  output[1](0).grad();

  EXPECT_NEAR(output[1](0).val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(a.adj(), -0.50107310888, 1e-5);
}

TEST(StanMathOde_ode_adams, std_vector_arg) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  std::vector<var> a = {1.5};

  var output = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a)[0][0];

  output.grad();

  EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(a[0].adj(), -0.50107310888, 1e-5);
}

TEST(StanMathOde_ode_adams, vector_arg) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  Eigen::Matrix<var, Eigen::Dynamic, 1> a(1);
  a << 1.5;

  var output = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a)[0][0];

  output.grad();

  EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(a(0).adj(), -0.50107310888, 1e-5);
}

TEST(StanMathOde_ode_adams, row_vector_arg) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  Eigen::Matrix<var, 1, Eigen::Dynamic> a(1);
  a << 1.5;

  var output = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a)[0][0];

  output.grad();

  EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(a(0).adj(), -0.50107310888, 1e-5);
}

TEST(StanMathOde_ode_adams, matrix_arg) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> a(1, 1);
  a << 1.5;

  var output = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a)[0][0];

  output.grad();

  EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(a(0, 0).adj(), -0.50107310888, 1e-5);
}

TEST(StanMathOde_ode_adams, scalar_std_vector_args) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  var a0 = 0.75;
  std::vector<var> a1 = {0.75};

  var output
      = stan::math::ode_adams(Cos2Arg(), y0, t0, ts, nullptr, a0, a1)[0][0];

  output.grad();

  EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(a0.adj(), -0.50107310888, 1e-5);
  EXPECT_NEAR(a1[0].adj(), -0.50107310888, 1e-5);
}

TEST(StanMathOde_ode_adams, std_vector_std_vector_args) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  var a0 = 1.5;
  std::vector<var> a1(1, a0);
  std::vector<std::vector<var>> a2(1, a1);

  var output = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a2)[0][0];

  output.grad();

  EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(a2[0][0].adj(), -0.50107310888, 1e-5);
}

TEST(StanMathOde_ode_adams, std_vector_vector_args) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  var a0 = 1.5;
  Eigen::Matrix<var, Eigen::Dynamic, 1> a1(1);
  a1 << a0;
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> a2(1, a1);

  var output = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a2)[0][0];

  output.grad();

  EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(a2[0](0).adj(), -0.50107310888, 1e-5);
}

TEST(StanMathOde_ode_adams, std_vector_row_vector_args) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  var a0 = 1.5;
  Eigen::Matrix<var, 1, Eigen::Dynamic> a1(1);
  a1 << a0;
  std::vector<Eigen::Matrix<var, 1, Eigen::Dynamic>> a2(1, a1);

  var output = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a2)[0][0];

  output.grad();

  EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(a2[0](0).adj(), -0.50107310888, 1e-5);
}

TEST(StanMathOde_ode_adams, std_vector_matrix_args) {
  using stan::math::var;

  Eigen::VectorXd y0 = Eigen::VectorXd::Zero(1);
  double t0 = 0.0;
  std::vector<double> ts = {1.1};

  var a0 = 1.5;
  Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic> a1(1, 1);
  a1 << a0;
  std::vector<Eigen::Matrix<var, Eigen::Dynamic, Eigen::Dynamic>> a2(1, a1);

  var output = stan::math::ode_adams(CosArg1(), y0, t0, ts, nullptr, a2)[0][0];

  output.grad();

  EXPECT_NEAR(output.val(), 0.66457668563, 1e-5);
  EXPECT_NEAR(a2[0](0).adj(), -0.50107310888, 1e-5);
}

struct ayt {
  template <typename T0, typename T_y, typename T2>
  inline Eigen::Matrix<stan::return_type_t<T_y, T2>, Eigen::Dynamic, 1>
  operator()(const T0& t, const T_y& y, std::ostream* msgs, const T2& a) const {
    return -a * y * t;
  }
};

TEST(StanMathOde_ode_adams, arg_combos_test) {
  var t0 = 0.5;
  var a = 0.2;
  std::vector<var> ts = {1.25};
  Eigen::Matrix<var, Eigen::Dynamic, 1> y0(1);
  y0 << 0.75;

  double t0d = stan::math::value_of(t0);
  double ad = stan::math::value_of(a);
  std::vector<double> tsd = stan::math::value_of(ts);
  Eigen::VectorXd y0d = stan::math::value_of(y0);

  auto check_yT = [&](auto yT) {
    EXPECT_NEAR(stan::math::value_of(yT),
                y0d(0) * exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)), 1e-5);
  };

  auto check_t0 = [&](var t0) {
    EXPECT_NEAR(
        t0.adj(),
        ad * t0d * y0d(0) * exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)),
        1e-5);
  };

  auto check_a = [&](var a) {
    EXPECT_NEAR(a.adj(),
                -0.5 * (tsd[0] * tsd[0] - t0d * t0d) * y0d(0)
                    * exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)),
                1e-5);
  };

  auto check_ts = [&](std::vector<var> ts) {
    EXPECT_NEAR(
        ts[0].adj(),
        -ad * tsd[0] * y0d(0) * exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)),
        1e-5);
  };

  auto check_y0 = [&](Eigen::Matrix<var, Eigen::Dynamic, 1> y0) {
    EXPECT_NEAR(y0(0).adj(), exp(-0.5 * ad * (tsd[0] * tsd[0] - t0d * t0d)),
                1e-5);
  };

  double yT1 = stan::math::ode_adams(ayt(), y0d, t0d, tsd, nullptr, ad)[0](0);
  check_yT(yT1);

  var yT2 = stan::math::ode_adams(ayt(), y0d, t0d, tsd, nullptr, a)[0](0);
  stan::math::set_zero_all_adjoints();
  yT2.grad();
  check_yT(yT2);
  check_a(a);

  var yT3 = stan::math::ode_adams(ayt(), y0d, t0d, ts, nullptr, ad)[0](0);
  stan::math::set_zero_all_adjoints();
  yT3.grad();
  check_yT(yT3);
  check_ts(ts);

  var yT4 = stan::math::ode_adams(ayt(), y0d, t0d, ts, nullptr, a)[0](0);
  stan::math::set_zero_all_adjoints();
  yT4.grad();
  check_yT(yT4);
  check_ts(ts);
  check_a(a);

  var yT5 = stan::math::ode_adams(ayt(), y0d, t0, tsd, nullptr, ad)[0](0);
  stan::math::set_zero_all_adjoints();
  yT5.grad();
  check_yT(yT5);
  check_t0(t0);

  var yT6 = stan::math::ode_adams(ayt(), y0d, t0, tsd, nullptr, a)[0](0);
  stan::math::set_zero_all_adjoints();
  yT6.grad();
  check_yT(yT6);
  check_t0(t0);
  check_a(a);

  var yT7 = stan::math::ode_adams(ayt(), y0d, t0, ts, nullptr, ad)[0](0);
  stan::math::set_zero_all_adjoints();
  yT7.grad();
  check_yT(yT7);
  check_t0(t0);
  check_ts(ts);

  var yT8 = stan::math::ode_adams(ayt(), y0d, t0, ts, nullptr, a)[0](0);
  stan::math::set_zero_all_adjoints();
  yT8.grad();
  check_yT(yT8);
  check_t0(t0);
  check_ts(ts);
  check_a(a);

  var yT9 = stan::math::ode_adams(ayt(), y0, t0d, tsd, nullptr, ad)[0](0);
  stan::math::set_zero_all_adjoints();
  yT9.grad();
  check_yT(yT9);
  check_y0(y0);

  var yT10 = stan::math::ode_adams(ayt(), y0, t0d, tsd, nullptr, a)[0](0);
  stan::math::set_zero_all_adjoints();
  yT10.grad();
  check_yT(yT10);
  check_y0(y0);
  check_a(a);

  var yT11 = stan::math::ode_adams(ayt(), y0, t0d, ts, nullptr, ad)[0](0);
  stan::math::set_zero_all_adjoints();
  yT11.grad();
  check_yT(yT11);
  check_y0(y0);
  check_ts(ts);

  var yT12 = stan::math::ode_adams(ayt(), y0, t0d, ts, nullptr, a)[0](0);
  stan::math::set_zero_all_adjoints();
  yT12.grad();
  check_yT(yT12);
  check_y0(y0);
  check_ts(ts);
  check_a(a);

  var yT13 = stan::math::ode_adams(ayt(), y0, t0, tsd, nullptr, ad)[0](0);
  stan::math::set_zero_all_adjoints();
  yT13.grad();
  check_yT(yT13);
  check_y0(y0);
  check_t0(t0);

  var yT14 = stan::math::ode_adams(ayt(), y0, t0, tsd, nullptr, a)[0](0);
  stan::math::set_zero_all_adjoints();
  yT14.grad();
  check_yT(yT14);
  check_y0(y0);
  check_t0(t0);
  check_a(a);

  var yT15 = stan::math::ode_adams(ayt(), y0, t0, ts, nullptr, ad)[0](0);
  stan::math::set_zero_all_adjoints();
  yT15.grad();
  check_yT(yT15);
  check_y0(y0);
  check_t0(t0);
  check_ts(ts);

  var yT16 = stan::math::ode_adams(ayt(), y0, t0, ts, nullptr, a)[0](0);
  stan::math::set_zero_all_adjoints();
  yT16.grad();
  check_yT(yT16);
  check_y0(y0);
  check_t0(t0);
  check_ts(ts);
  check_a(a);
}
