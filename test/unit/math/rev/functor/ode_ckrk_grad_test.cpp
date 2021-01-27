#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <stan/math/rev/functor/gradient.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>

using stan::math::ode_ckrk;
using std::cos;
using std::sin;

double y1(double t, double omega, double chi) { return chi * cos(omega * t); }

double dy1_domega(double t, double omega, double chi) {
  return -t * chi * sin(omega * t);
}

double dy1_dchi(double t, double omega, double chi) { return cos(omega * t); }

double y2(double t, double omega, double chi) {
  return -omega * chi * sin(omega * t);
}

double dy2_domega(double t, double omega, double chi) {
  return -chi * (sin(omega * t) + omega * t * cos(omega * t));
}

double dy2_dchi(double t, double omega, double chi) {
  return -omega * sin(omega * t);
}

class sho_functor {
 public:
  template <typename T0, typename T1, typename T2>
  inline Eigen::Matrix<stan::return_type_t<T1, T2>, -1, 1> operator()(
      const T0& t_in, const Eigen::Matrix<T1, -1, 1>& y_in, std::ostream* msgs,
      const std::vector<T2>& theta, const std::vector<double>& x,
      const std::vector<int>& x_int) const {
    if (y_in.size() != 2)
      throw std::domain_error("Functor called with inconsistent state");

    Eigen::Matrix<stan::return_type_t<T1, T2>, -1, 1> f(2);
    f << (y_in(1)), (-theta.at(0) * theta.at(0) * y_in(0));

    return f;
  }
};

class test_functor_double_var_1 {
 public:
  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    std::vector<T> theta;
    theta.push_back(x(0));

    Eigen::Matrix<double, -1, 1> y0(2);
    y0 << 1.25, 0.0;

    double t0 = 0.0;
    std::vector<double> ts;
    ts.push_back(5.0);

    std::vector<double> data;
    std::vector<int> data_int;

    auto ys = ode_ckrk(sho, y0, t0, ts, 0, theta, data, data_int);

    return ys[0][0];
  }
};

class test_functor_double_var_2 {
 public:
  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    std::vector<T> theta;
    theta.push_back(x(0));

    Eigen::Matrix<double, -1, 1> y0(2);
    y0 << 1.25, 0.0;

    double t0 = 0.0;
    std::vector<double> ts;
    ts.push_back(5.0);

    std::vector<double> data;
    std::vector<int> data_int;

    auto ys = ode_ckrk(sho, y0, t0, ts, 0, theta, data, data_int);

    return ys[0][1];
  }
};

TEST(StanMathOdeIntegrateODEGradArr, ckrk_double_var) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;

  Eigen::VectorXd x(1);
  x(0) = omega;

  double f;
  Eigen::VectorXd grad(1);

  test_functor_double_var_1 func1;
  stan::math::gradient(func1, x, f, grad);

  EXPECT_NEAR(y1(t, omega, chi), f, 1e-5);
  EXPECT_NEAR(dy1_domega(t, omega, chi), grad(0), 1e-5);

  test_functor_double_var_2 func2;
  stan::math::gradient(func2, x, f, grad);

  EXPECT_NEAR(y2(t, omega, chi), f, 1e-5);
  EXPECT_NEAR(dy2_domega(t, omega, chi), grad(0), 1e-5);
}

class test_functor_var_double_1 {
 public:
  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    std::vector<double> theta;
    theta.push_back(0.5);

    Eigen::Matrix<T, -1, 1> y0(2);
    y0 << x(0), 0.0;

    double t0 = 0.0;
    std::vector<double> ts;
    ts.push_back(5.0);

    std::vector<double> data;
    std::vector<int> data_int;

    auto ys = ode_ckrk(sho, y0, t0, ts, 0, theta, data, data_int);

    return ys[0][0];
  }
};

class test_functor_var_double_2 {
 public:
  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    std::vector<double> theta;
    theta.push_back(0.5);

    Eigen::Matrix<T, -1, 1> y0(2);
    y0 << x(0), 0.0;

    double t0 = 0.0;
    std::vector<double> ts;
    ts.push_back(5.0);

    std::vector<double> data;
    std::vector<int> data_int;

    auto ys = ode_ckrk(sho, y0, t0, ts, 0, theta, data, data_int);

    return ys[0][1];
  }
};

TEST(StanMathOdeIntegrateODEGradArr, ckrk_var_double) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;

  Eigen::VectorXd x(1);
  x(0) = chi;

  double f;
  Eigen::VectorXd grad(1);

  test_functor_var_double_1 func1;
  stan::math::gradient(func1, x, f, grad);

  EXPECT_NEAR(y1(t, omega, chi), f, 1e-5);
  EXPECT_NEAR(dy1_dchi(t, omega, chi), grad(0), 1e-5);

  test_functor_var_double_2 func2;
  stan::math::gradient(func2, x, f, grad);

  EXPECT_NEAR(y2(t, omega, chi), f, 1e-5);
  EXPECT_NEAR(dy2_dchi(t, omega, chi), grad(0), 1e-5);
}

class test_functor_var_var_1 {
 public:
  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    std::vector<T> theta;
    theta.push_back(x(0));

    Eigen::Matrix<T, -1, 1> y0(2);
    y0 << x(1), 0.0;

    double t0 = 0.0;
    std::vector<double> ts;
    ts.push_back(5.0);

    std::vector<double> data;
    std::vector<int> data_int;

    auto ys = ode_ckrk(sho, y0, t0, ts, 0, theta, data, data_int);

    return ys[0][0];
  }
};

class test_functor_var_var_2 {
 public:
  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    std::vector<T> theta;
    theta.push_back(x(0));

    Eigen::Matrix<T, -1, 1> y0(2);
    y0 << x(1), 0.0;

    double t0 = 0.0;
    std::vector<double> ts;
    ts.push_back(5.0);

    std::vector<double> data;
    std::vector<int> data_int;

    auto ys = ode_ckrk(sho, y0, t0, ts, 0, theta, data, data_int);

    return ys[0][1];
  }
};

TEST(StanMathOdeIntegrateODEGradArr, ckrk_var_var) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;

  Eigen::VectorXd x(2);
  x(0) = omega;
  x(1) = chi;

  double f;
  Eigen::VectorXd grad(2);

  test_functor_var_var_1 func1;
  stan::math::gradient(func1, x, f, grad);

  EXPECT_NEAR(y1(t, omega, chi), f, 1e-5);
  EXPECT_NEAR(dy1_domega(t, omega, chi), grad(0), 1e-5);
  EXPECT_NEAR(dy1_dchi(t, omega, chi), grad(1), 1e-5);

  test_functor_var_var_2 func2;
  stan::math::gradient(func2, x, f, grad);

  EXPECT_NEAR(y2(t, omega, chi), f, 1e-5);
  EXPECT_NEAR(dy2_domega(t, omega, chi), grad(0), 1e-5);
  EXPECT_NEAR(dy2_dchi(t, omega, chi), grad(1), 1e-5);
}
