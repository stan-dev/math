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

struct sho_functor {
  template <typename T0, typename T1, typename T2>
  inline Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1>
  operator()(const T0& t, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             std::ostream* msgs, const T2& omega) const {
    Eigen::Matrix<stan::return_type_t<T1, T2>, Eigen::Dynamic, 1> out(2);
    out << y.coeff(1), -omega * omega * y.coeff(0);
    return out;
  }
};

template <int state>
class test_functor_double_var {
  int lmm_;

 public:
  explicit test_functor_double_var(int lmm) : lmm_(lmm) {}

  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    T omega = x(0);

    Eigen::VectorXd y0(2);
    y0 << 1.25, 0.0;

    double t0 = 0.0;
    std::vector<double> ts{2.5, 5.0};

    std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> ys
        = (lmm_ == TEST_CVODES_ADAMS)
              ? stan::math::ode_bdf_adjoint(sho, y0, t0, ts, nullptr, omega)
              : stan::math::ode_bdf(sho, y0, t0, ts, nullptr, omega);

    return ys[1](state);
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

  stan::math::nested_rev_autodiff nested;

  {  // Adams
    test_functor_double_var<0> func1(TEST_CVODES_ADAMS);
    stan::math::gradient(func1, x, f, grad);

    EXPECT_NEAR(y1(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy1_domega(t, omega, chi), grad(0), 1e-6);

    test_functor_double_var<1> func2(TEST_CVODES_ADAMS);
    stan::math::gradient(func2, x, f, grad);

    EXPECT_NEAR(y2(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy2_domega(t, omega, chi), grad(0), 1e-6);
  }

  {  // bdf
    test_functor_double_var<0> func1(TEST_CVODES_BDF);
    stan::math::gradient(func1, x, f, grad);

    EXPECT_NEAR(y1(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy1_domega(t, omega, chi), grad(0), 1e-7);

    test_functor_double_var<1> func2(TEST_CVODES_BDF);
    stan::math::gradient(func2, x, f, grad);

    EXPECT_NEAR(y2(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy2_domega(t, omega, chi), grad(0), 1e-7);
  }
}

template <int state>
class test_functor_var_double {
  int lmm_;

 public:
  explicit test_functor_var_double(int lmm) : lmm_(lmm) {}

  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    double omega = 0.5;

    Eigen::Matrix<T, Eigen::Dynamic, 1> y0(2);
    y0 << x(0), 0.0;

    double t0 = 0.0;
    std::vector<double> ts{5.0};

    std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> ys
        = (lmm_ == TEST_CVODES_ADAMS)
              ? stan::math::ode_bdf_adjoint(sho, y0, t0, ts, nullptr, omega)
              : stan::math::ode_bdf(sho, y0, t0, ts, nullptr, omega);

    return ys[0](state);
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

  stan::math::nested_rev_autodiff nested;
  {  // adams
    test_functor_var_double<0> func1(TEST_CVODES_ADAMS);
    stan::math::gradient(func1, x, f, grad);

    EXPECT_NEAR(y1(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy1_dchi(t, omega, chi), grad(0), 1e-7);

    test_functor_var_double<1> func2(TEST_CVODES_ADAMS);
    stan::math::gradient(func2, x, f, grad);

    EXPECT_NEAR(y2(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy2_dchi(t, omega, chi), grad(0), 1e-7);
  }

  {  // bdf
    test_functor_var_double<0> func1(TEST_CVODES_BDF);
    stan::math::gradient(func1, x, f, grad);

    EXPECT_NEAR(y1(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy1_dchi(t, omega, chi), grad(0), 1e-7);

    test_functor_var_double<1> func2(TEST_CVODES_BDF);
    stan::math::gradient(func2, x, f, grad);

    EXPECT_NEAR(y2(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy2_dchi(t, omega, chi), grad(0), 1e-7);
  }
}

template <int state>
class test_functor_var_var {
  int lmm_;

 public:
  explicit test_functor_var_var(int lmm) : lmm_(lmm) {}

  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    T omega = x(0);

    Eigen::Matrix<T, Eigen::Dynamic, 1> y0(2);
    y0 << x(1), 0.0;

    double t0 = 0.0;
    std::vector<double> ts{2.0, 5.0};

    std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> ys
        = (lmm_ == TEST_CVODES_ADAMS)
              ? stan::math::ode_adams_tol(sho, y0, t0, ts, 1E-10, 1E-10, 10000,
                                          nullptr, omega)
              : stan::math::ode_bdf(sho, y0, t0, ts, nullptr, omega);

    return ys[1](state);
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
    stan::math::nested_rev_autodiff nested;

    test_functor_var_var<0> func1(TEST_CVODES_ADAMS);
    stan::math::gradient(func1, x, f, grad);

    EXPECT_NEAR(y1(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy1_domega(t, omega, chi), grad(0), 1e-7);
    EXPECT_NEAR(dy1_dchi(t, omega, chi), grad(1), 1e-7);

    test_functor_var_var<1> func2(TEST_CVODES_ADAMS);
    stan::math::gradient(func2, x, f, grad);

    EXPECT_NEAR(y2(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy2_domega(t, omega, chi), grad(0), 1e-6);
    EXPECT_NEAR(dy2_dchi(t, omega, chi), grad(1), 1e-7);
  }

  {
    stan::math::nested_rev_autodiff nested;

    test_functor_var_var<0> func1(TEST_CVODES_BDF);
    stan::math::gradient(func1, x, f, grad);

    EXPECT_NEAR(y1(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy1_domega(t, omega, chi), grad(0), 1e-7);
    EXPECT_NEAR(dy1_dchi(t, omega, chi), grad(1), 1e-7);

    test_functor_var_var<1> func2(TEST_CVODES_BDF);
    stan::math::gradient(func2, x, f, grad);

    EXPECT_NEAR(y2(t, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy2_domega(t, omega, chi), grad(0), 1e-6);
    EXPECT_NEAR(dy2_dchi(t, omega, chi), grad(1), 1e-7);
  }
}



template <int state>
class test_functor_sum_var_var {
  int lmm_;

 public:
  explicit test_functor_sum_var_var(int lmm) : lmm_(lmm) {}

  template <typename T>
  inline T operator()(Eigen::Matrix<T, Eigen::Dynamic, 1>& x) const {
    sho_functor sho;

    T omega = x(0);

    Eigen::Matrix<T, Eigen::Dynamic, 1> y0(2);
    y0 << x(1), 0.0;

    double t0 = 0.0;
    std::vector<double> ts{2.0, 5.0};

    std::vector<Eigen::Matrix<T, Eigen::Dynamic, 1>> ys
        = (lmm_ == TEST_CVODES_ADAMS)
              ? stan::math::ode_adams_tol(sho, y0, t0, ts, 1E-10, 1E-10, 10000,
                                          nullptr, omega)
              : stan::math::ode_bdf(sho, y0, t0, ts, nullptr, omega);

    return stan::math::sum(ys[0](state) + ys[1](state));
  }
};

TEST(StanMathOdeIntegrateODEGradMat, sum_var_var) {
  double omega = 0.5;
  double chi = 1.25;
  double t1 = 2.0;
  double t2 = 5;

  Eigen::VectorXd x(2);
  x(0) = omega;
  x(1) = chi;

  double f;
  Eigen::VectorXd grad(2);

  {
    stan::math::nested_rev_autodiff nested;

    test_functor_sum_var_var<0> func1(TEST_CVODES_ADAMS);
    stan::math::gradient(func1, x, f, grad);

    EXPECT_NEAR(y1(t1, omega, chi) + y1(t2, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy1_domega(t1, omega, chi) + dy1_domega(t2, omega, chi), grad(0), 1e-7);
    EXPECT_NEAR(dy1_dchi(t1, omega, chi) + dy1_dchi(t2, omega, chi), grad(1), 1e-7);

    test_functor_sum_var_var<1> func2(TEST_CVODES_ADAMS);
    stan::math::gradient(func2, x, f, grad);

    EXPECT_NEAR(y2(t1, omega, chi) + y2(t2, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy2_domega(t1, omega, chi) + dy2_domega(t2, omega, chi), grad(0), 1e-6);
    EXPECT_NEAR(dy2_dchi(t1, omega, chi) + dy2_dchi(t2, omega, chi), grad(1), 1e-7);
  }

  {
    stan::math::nested_rev_autodiff nested;

    test_functor_sum_var_var<0> func1(TEST_CVODES_BDF);
    stan::math::gradient(func1, x, f, grad);

    EXPECT_NEAR(y1(t1, omega, chi) + y1(t2, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy1_domega(t1, omega, chi) + dy1_domega(t2, omega, chi), grad(0), 1e-7);
    EXPECT_NEAR(dy1_dchi(t1, omega, chi) + dy1_dchi(t2, omega, chi), grad(1), 1e-7);

    test_functor_sum_var_var<1> func2(TEST_CVODES_BDF);
    stan::math::gradient(func2, x, f, grad);

    EXPECT_NEAR(y2(t1, omega, chi) + y2(t2, omega, chi), f, 1e-8);
    EXPECT_NEAR(dy2_domega(t1, omega, chi) + dy2_domega(t2, omega, chi), grad(0), 1e-6);
    EXPECT_NEAR(dy2_dchi(t1, omega, chi) + dy2_dchi(t2, omega, chi), grad(1), 1e-7);
  }
}
