#ifndef STAN_MATH_TEST_SHO_GRAD_TEST_FIXTURE_HPP
#define STAN_MATH_TEST_SHO_GRAD_TEST_FIXTURE_HPP

// #include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>

struct sho_grad_ode_test : public testing::Test {
  class sho_functor {
   public:
    template <typename T0, typename T1, typename T2>
    inline Eigen::Matrix<stan::return_type_t<T1, T2>, -1, 1> operator()(
        const T0& t_in, const Eigen::Matrix<T1, -1, 1>& y_in,
        std::ostream* msgs, const std::vector<T2>& theta,
        const std::vector<double>& x, const std::vector<int>& x_int) const {
      if (y_in.size() != 2)
        throw std::domain_error("Functor called with inconsistent state");

      Eigen::Matrix<stan::return_type_t<T1, T2>, -1, 1> f(2);
      f << (y_in(1)), (-theta.at(0) * theta.at(0) * y_in(0));

      return f;
    }
  };

  template <typename solver_functor_t>
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

      solver_functor_t sol;
      auto ys = sol(sho, y0, t0, ts, 0, theta, data, data_int);

      return ys[0][0];
    }
  };

  template <typename solver_functor_t>
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

      solver_functor_t sol;
      auto ys = sol(sho, y0, t0, ts, 0, theta, data, data_int);

      return ys[0][1];
    }
  };

  template <typename solver_functor_t>
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

      solver_functor_t sol;
      auto ys = sol(sho, y0, t0, ts, 0, theta, data, data_int);

      return ys[0][0];
    }
  };

  template <typename solver_functor_t>
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

      solver_functor_t sol;
      auto ys = sol(sho, y0, t0, ts, 0, theta, data, data_int);

      return ys[0][1];
    }
  };

  template <typename solver_functor_t>
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

      solver_functor_t sol;
      auto ys = sol(sho, y0, t0, ts, 0, theta, data, data_int);

      return ys[0][0];
    }
  };

  template <typename solver_functor_t>
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

      solver_functor_t sol;
      auto ys = sol(sho, y0, t0, ts, 0, theta, data, data_int);

      return ys[0][1];
    }
  };

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

  void SetUp() {
    // make sure memory's clean before starting each test
    stan::math::recover_memory();
  }

  sho_grad_ode_test() { SetUp(); }
};

#endif
