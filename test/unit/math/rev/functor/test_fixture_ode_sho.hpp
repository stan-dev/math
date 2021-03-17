#include <stan/math/rev.hpp>
#include <boost/numeric/odeint.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/functor/harmonic_oscillator.hpp>
#include <test/unit/math/rev/functor/coupled_mm.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

struct harmonic_oscillator_ode_base {
  harm_osc_ode_fun f;
  harm_osc_ode_fun_eigen f_eigen;
  harm_osc_ode_data_fun f_data;
  harm_osc_ode_data_fun_eigen f_data_eigen;
  harm_osc_ode_wrong_size_1_fun f_wrong_size_1;
  harm_osc_ode_wrong_size_2_fun f_wrong_size_2;

  std::vector<double> theta;
  Eigen::VectorXd y0;
  std::vector<double> y0_eigen;
  double t0;
  std::vector<double> ts;
  std::vector<double> x_r;
  std::vector<int> x_i;
  double rtol;
  double atol;
  int max_num_step;

  harmonic_oscillator_ode_base() :
    theta{0.15}, y0(2), t0(0), ts(100),
    rtol(1.e-6), atol(1.e-8), max_num_step(100000)
  {
    y0 << 1.0, 0.0;
    for (size_t i = 0; i < ts.size(); ++i) {
      ts[i] = t0 + 0.1 * (i + 1);
    }
  }
};

template<typename T>
struct harmonic_oscillator_test : public harmonic_oscillator_ode_base,
                                         public ODETestFixture<harmonic_oscillator_test<T>> {
  harmonic_oscillator_test() : harmonic_oscillator_ode_base()
  {}

  auto apply_solver() {
    std::tuple_element_t<0, T> sol;
    return sol(f_eigen, y0, t0, ts, nullptr, theta, x_r, x_i);
  }

  auto apply_solver_tol() {
    std::tuple_element_t<1, T> sol;
    sol(f_eigen, y0, t0, ts, rtol, atol, max_num_step, nullptr, theta, x_r, x_i);
  }

  void test_good() {
    ASSERT_NO_THROW(apply_solver());
  }

  void test_bad() {
    const Eigen::VectorXd y0_(y0);

    y0 = Eigen::VectorXd();
    EXPECT_THROW_MSG(apply_solver(), std::invalid_argument, "initial state has size 0");
    y0 = y0_;

    const double t0_ = t0;
    t0 = 2.0;
    EXPECT_THROW_MSG(apply_solver(), std::domain_error, "initial time is 2, but must be less than 0.1");
    t0 = t0_;

    const std::vector<double> ts_ = ts;
    ts = std::vector<double>();
    EXPECT_THROW_MSG(apply_solver(), std::invalid_argument, "times has size 0");
    ts = ts_;

    ts = std::vector<double>{3, 1};
    EXPECT_THROW_MSG(apply_solver(), std::domain_error, "times is not a valid sorted vector");
    ts = ts_;

    const std::vector<double> theta_= theta;
    const std::vector<double> x_r_ = x_r;
    const std::vector<int> x_i_ = x_i;

    const double rtol_ = rtol;
    rtol = -1;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "relative_tolerance");
    rtol = rtol_;

    const double atol_ = atol;
    atol = -1;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "absolute_tolerance");
    atol = atol_;

    // NaN errors
    double nan = std::numeric_limits<double>::quiet_NaN();
    std::stringstream expected_is_nan;
    expected_is_nan << "is " << nan;

    y0[0] = nan;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial state");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_nan.str());
    y0 = y0_;

    t0 = nan;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial time");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_nan.str());
    t0 = t0_;

    ts[0] = nan;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "times");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_nan.str());
    ts = ts_;

    theta[0] = nan;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_nan.str());
    theta = theta_;

    x_r.push_back(nan);
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_nan.str());
    x_r = x_r_;

    // inf test
    std::stringstream expected_is_inf;
    expected_is_inf << "is " << std::numeric_limits<double>::infinity();
    std::stringstream expected_is_neg_inf;
    expected_is_neg_inf << "is " << -std::numeric_limits<double>::infinity();    
    double inf = std::numeric_limits<double>::infinity();

    y0[0] = inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial state");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_inf.str());    
    y0[0] = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial state");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_neg_inf.str());    
    y0 = y0_;
    
    t0 = inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial time");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_inf.str());    
    t0 = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial time");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_neg_inf.str());    
    t0 = t0_;

    ts.back() = inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "times");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_inf.str());    
    ts.back() = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "times");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_neg_inf.str());    
    ts = ts_;

    theta[0] = inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_inf.str());
    theta[0] = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_neg_inf.str());
    theta = theta_;

    x_r = std::vector<double>{inf};
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_inf.str());
    x_r[0] = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "ode parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, expected_is_neg_inf.str());
    x_r = x_r_;
  }
};

template<typename T>
struct harmonic_oscillator_bad_ode_test : public harmonic_oscillator_ode_base,
                                           public ODETestFixture<harmonic_oscillator_bad_ode_test<T>> {
  harmonic_oscillator_bad_ode_test() : harmonic_oscillator_ode_base()
  {}

  auto apply_solver() {
    std::tuple_element_t<0, T> sol;
    return sol(f_wrong_size_1, stan::math::to_array_1d(y0), t0, ts, theta, x_r, x_i, 0);
  }

  auto apply_solver_tol() {
    std::tuple_element_t<1, T> sol;
    sol(f_wrong_size_1, stan::math::to_array_1d(y0), t0, ts, theta, x_r, x_i, 0, rtol, atol, max_num_step);
  }

  void test_bad_ode() {
    std::string error_msg = "dy_dt (3) and states (2) must match in size";
    EXPECT_THROW_MSG(apply_solver_tol(), std::invalid_argument, error_msg);
  }
};

template<typename T>
struct harmonic_oscillator_data_test : public harmonic_oscillator_ode_base,
                                            public ODETestFixture<harmonic_oscillator_data_test<T>> {
  harmonic_oscillator_data_test() : harmonic_oscillator_ode_base()
  {}

  auto apply_solver() {
    std::tuple_element_t<0, T> sol;
    return sol(f_data, stan::math::to_array_1d(y0), t0, ts, theta, x_r, x_i, 0);
  }

  auto apply_solver_tol() {
    std::tuple_element_t<1, T> sol;
    sol(f_data, stan::math::to_array_1d(y0), t0, ts, theta, x_r, x_i, 0, rtol, atol, max_num_step);
  }

  void test_bad_param_and_data() {
    const std::vector<double> theta_= theta;
    theta = std::vector<double>();
    EXPECT_THROW_MSG(apply_solver(), std::out_of_range, "vector");
    theta = theta_;

    const std::vector<double> x_r_ = x_r;
    x_r = std::vector<double>();
    EXPECT_THROW_MSG(apply_solver(), std::out_of_range, "vector");
    x_r = x_r_;

    const std::vector<int> x_i_ = x_i;
    x_i = std::vector<int>();
    EXPECT_THROW_MSG(apply_solver(), std::out_of_range, "vector");

    max_num_step = -1;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "max_num_steps");
  }
};
