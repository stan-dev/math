#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <stan/math/rev/functor/gradient.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>

#define TEST_CVODES_ADAMS 1
#define TEST_CVODES_BDF 2

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
  inline std::vector<stan::return_type_t<T1, T2>> operator()(
      const T0& t_in,                 // time
      const std::vector<T1>& y_in,    // state
      const std::vector<T2>& theta,   // parameters
      const std::vector<double>& x,   // double data
      const std::vector<int>& x_int,  // integer data
      std::ostream* msgs) const {     // error stream
    if (y_in.size() != 2)
      throw std::domain_error("Functor called with inconsistent state");

    std::vector<stan::return_type_t<T1, T2>> f;
    f.push_back(y_in.at(1));
    f.push_back(-theta.at(0) * theta.at(0) * y_in.at(0));

    return f;
  }
};

class test_functor_double_var_1 {
  int lmm_;

 public:
  explicit test_functor_double_var_1(int lmm) : lmm_(lmm) {}

  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    std::vector<T> theta;
    theta.push_back(x(0));

    std::vector<double> y0;
    y0.push_back(1.25);
    y0.push_back(0.0);

    double t0 = 0.0;
    std::vector<double> ts;
    ts.push_back(5.0);

    std::vector<double> data;
    std::vector<int> data_int;

    std::vector<std::vector<T>> ys
        = (lmm_ == TEST_CVODES_ADAMS
               ? stan::math::integrate_ode_adams(sho, y0, t0, ts, theta, data,
                                                 data_int)
               : stan::math::integrate_ode_bdf(sho, y0, t0, ts, theta, data,
                                               data_int));

    return ys[0][0];
  }
};

class test_functor_double_var_2 {
  int lmm_;

 public:
  explicit test_functor_double_var_2(int lmm) : lmm_(lmm) {}

  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    std::vector<T> theta;
    theta.push_back(x(0));

    std::vector<double> y0;
    y0.push_back(1.25);
    y0.push_back(0.0);

    double t0 = 0.0;
    std::vector<double> ts;
    ts.push_back(5.0);

    std::vector<double> data;
    std::vector<int> data_int;

    std::vector<std::vector<T>> ys
        = (lmm_ == TEST_CVODES_ADAMS
               ? stan::math::integrate_ode_adams(sho, y0, t0, ts, theta, data,
                                                 data_int)
               : stan::math::integrate_ode_bdf(sho, y0, t0, ts, theta, data,
                                               data_int));

    return ys[0][1];
  }
};

TEST(StanMathOdeIntegrateODEGradMat, double_var) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;

  Eigen::VectorXd x(1);
  x(0) = omega;

  double f;
  Eigen::VectorXd grad(1);

  {  // Adams
    test_functor_double_var_1 func1(TEST_CVODES_ADAMS);
    stan::math::gradient(func1, x, f, grad);

    EXPECT_NEAR(y1(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy1_domega(t, omega, chi), grad(0), 1e-7);

    test_functor_double_var_2 func2(TEST_CVODES_ADAMS);
    stan::math::gradient(func2, x, f, grad);

    EXPECT_NEAR(y2(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy2_domega(t, omega, chi), grad(0), 1e-7);
  }

  {  // bdf
    test_functor_double_var_1 func1(TEST_CVODES_BDF);
    stan::math::gradient(func1, x, f, grad);

    EXPECT_NEAR(y1(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy1_domega(t, omega, chi), grad(0), 1e-7);

    test_functor_double_var_2 func2(TEST_CVODES_BDF);
    stan::math::gradient(func2, x, f, grad);

    EXPECT_NEAR(y2(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy2_domega(t, omega, chi), grad(0), 1e-7);
  }
}

class test_functor_var_double_1 {
  int lmm_;

 public:
  explicit test_functor_var_double_1(int lmm) : lmm_(lmm) {}

  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    std::vector<double> theta;
    theta.push_back(0.5);

    std::vector<T> y0;
    y0.push_back(x(0));
    y0.push_back(0.0);

    double t0 = 0.0;
    std::vector<double> ts;
    ts.push_back(5.0);

    std::vector<double> data;
    std::vector<int> data_int;

    std::vector<std::vector<T>> ys
        = (lmm_ == TEST_CVODES_ADAMS
               ? stan::math::integrate_ode_adams(sho, y0, t0, ts, theta, data,
                                                 data_int)
               : stan::math::integrate_ode_bdf(sho, y0, t0, ts, theta, data,
                                               data_int));

    return ys[0][0];
  }
};

class test_functor_var_double_2 {
  int lmm_;

 public:
  explicit test_functor_var_double_2(int lmm) : lmm_(lmm) {}

  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    std::vector<double> theta;
    theta.push_back(0.5);

    std::vector<T> y0;
    y0.push_back(x(0));
    y0.push_back(0.0);

    double t0 = 0.0;
    std::vector<double> ts;
    ts.push_back(5.0);

    std::vector<double> data;
    std::vector<int> data_int;

    std::vector<std::vector<T>> ys
        = (lmm_ == TEST_CVODES_ADAMS
               ? stan::math::integrate_ode_adams(sho, y0, t0, ts, theta, data,
                                                 data_int)
               : stan::math::integrate_ode_bdf(sho, y0, t0, ts, theta, data,
                                               data_int));

    return ys[0][1];
  }
};

TEST(StanMathOdeIntegrateODEGradMat, var_double) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;

  Eigen::VectorXd x(1);
  x(0) = chi;

  double f;
  Eigen::VectorXd grad(1);

  {  // adams
    test_functor_var_double_1 func1(TEST_CVODES_ADAMS);
    stan::math::gradient(func1, x, f, grad);

    EXPECT_NEAR(y1(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy1_dchi(t, omega, chi), grad(0), 1e-7);

    test_functor_var_double_2 func2(TEST_CVODES_ADAMS);
    stan::math::gradient(func2, x, f, grad);

    EXPECT_NEAR(y2(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy2_dchi(t, omega, chi), grad(0), 1e-7);
  }

  {  // bdf
    test_functor_var_double_1 func1(TEST_CVODES_BDF);
    stan::math::gradient(func1, x, f, grad);

    EXPECT_NEAR(y1(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy1_dchi(t, omega, chi), grad(0), 1e-7);

    test_functor_var_double_2 func2(TEST_CVODES_BDF);
    stan::math::gradient(func2, x, f, grad);

    EXPECT_NEAR(y2(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy2_dchi(t, omega, chi), grad(0), 1e-7);
  }
}

class test_functor_var_var_1 {
  int lmm_;

 public:
  explicit test_functor_var_var_1(int lmm) : lmm_(lmm) {}

  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    std::vector<T> theta;
    theta.push_back(x(0));

    std::vector<T> y0;
    y0.push_back(x(1));
    y0.push_back(0.0);

    double t0 = 0.0;
    std::vector<double> ts;
    ts.push_back(5.0);

    std::vector<double> data;
    std::vector<int> data_int;

    std::vector<std::vector<T>> ys
        = (lmm_ == TEST_CVODES_ADAMS
               ? stan::math::integrate_ode_adams(sho, y0, t0, ts, theta, data,
                                                 data_int)
               : stan::math::integrate_ode_bdf(sho, y0, t0, ts, theta, data,
                                               data_int));

    return ys[0][0];
  }
};

class test_functor_var_var_2 {
  int lmm_;

 public:
  explicit test_functor_var_var_2(int lmm) : lmm_(lmm) {}

  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    std::vector<T> theta;
    theta.push_back(x(0));

    std::vector<T> y0;
    y0.push_back(x(1));
    y0.push_back(0.0);

    double t0 = 0.0;
    std::vector<double> ts;
    ts.push_back(5.0);

    std::vector<double> data;
    std::vector<int> data_int;

    std::vector<std::vector<T>> ys
        = (lmm_ == TEST_CVODES_ADAMS
               ? stan::math::integrate_ode_adams(sho, y0, t0, ts, theta, data,
                                                 data_int)
               : stan::math::integrate_ode_bdf(sho, y0, t0, ts, theta, data,
                                               data_int));

    return ys[0][1];
  }
};

TEST(StanMathOdeIntegrateODEGradMat, var_var) {
  double omega = 0.5;
  double chi = 1.25;
  double t = 5;

  Eigen::VectorXd x(2);
  x(0) = omega;
  x(1) = chi;

  double f;
  Eigen::VectorXd grad(2);

  {
    test_functor_var_var_1 func1(TEST_CVODES_ADAMS);
    stan::math::gradient(func1, x, f, grad);

    EXPECT_NEAR(y1(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy1_domega(t, omega, chi), grad(0), 1e-7);
    EXPECT_NEAR(dy1_dchi(t, omega, chi), grad(1), 1e-7);

    test_functor_var_var_2 func2(TEST_CVODES_ADAMS);
    stan::math::gradient(func2, x, f, grad);

    EXPECT_NEAR(y2(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy2_domega(t, omega, chi), grad(0), 1e-7);
    EXPECT_NEAR(dy2_dchi(t, omega, chi), grad(1), 1e-7);
  }

  {
    test_functor_var_var_1 func1(TEST_CVODES_BDF);
    stan::math::gradient(func1, x, f, grad);

    EXPECT_NEAR(y1(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy1_domega(t, omega, chi), grad(0), 1e-7);
    EXPECT_NEAR(dy1_dchi(t, omega, chi), grad(1), 1e-7);

    test_functor_var_var_2 func2(TEST_CVODES_BDF);
    stan::math::gradient(func2, x, f, grad);

    EXPECT_NEAR(y2(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy2_domega(t, omega, chi), grad(0), 1e-7);
    EXPECT_NEAR(dy2_dchi(t, omega, chi), grad(1), 1e-7);
  }
}
