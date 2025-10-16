#ifndef STAN_MATH_TEST_FIXTURE_DAE_CHEM_HPP
#define STAN_MATH_TEST_FIXTURE_DAE_CHEM_HPP

#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/functor/coupled_mm.hpp>
#include <test/unit/math/rev/functor/test_fixture_ode.hpp>
#include <test/unit/util.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <string>

template <typename T>
struct chemical_kinetics_dae_base {
  struct chemical_kinetics_fun {
    template <typename T0, typename Tyy, typename Typ, typename Tpar>
    inline Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> operator()(
        const T0& t_in, const Eigen::Matrix<Tyy, -1, 1>& yy,
        const Eigen::Matrix<Typ, -1, 1>& yp, std::ostream* msgs,
        const std::vector<Tpar>& theta, const std::vector<double>& x_r,
        const std::vector<int>& x_i) const {
      if (yy.size() != 3 || yp.size() != 3)
        throw std::domain_error(
            "this function was called with inconsistent state");

      Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> res(3);

      auto yy1 = yy(0);
      auto yy2 = yy(1);
      auto yy3 = yy(2);

      auto yp1 = yp(0);
      auto yp2 = yp(1);
      // auto yp3 = yp.at(2);

      auto p1 = theta.at(0);
      auto p2 = theta.at(1);
      auto p3 = theta.at(2);

      res[0] = yp1 + p1 * yy1 - p2 * yy2 * yy3;
      res[1] = yp2 - p1 * yy1 + p2 * yy2 * yy3 + p3 * yy2 * yy2;
      res[2] = yy1 + yy2 + yy3 - 1.0;

      return res;
    }
  };

  struct chemical_kinetics_data_fun {
    template <typename T0, typename Tyy, typename Typ, typename Tpar>
    inline Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> operator()(
        const T0& t_in, const Eigen::Matrix<Tyy, -1, 1>& yy,
        const Eigen::Matrix<Typ, -1, 1>& yp, std::ostream* msgs,
        const std::vector<Tpar>& theta, const std::vector<double>& x_r,
        const std::vector<int>& x_i) const {
      if (yy.size() != 3 || yp.size() != 3)
        throw std::domain_error(
            "this function was called with inconsistent state");

      Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> res(3);

      auto yy1 = yy(0);
      auto yy2 = yy(1);
      auto yy3 = yy(2);

      auto yp1 = yp(0);
      auto yp2 = yp(1);
      // auto yp3 = yp.at(2);

      auto p1 = theta.at(0);
      auto p2 = x_r.at(0);
      auto p3 = x_r.at(1);

      res[0] = yp1 + p1 * yy1 - p2 * yy2 * yy3;
      res[1] = yp2 - p1 * yy1 + p2 * yy2 * yy3 + p3 * yy2 * yy2;
      res[2] = yy1 + yy2 + yy3 - 1.0;

      return res;
    }
  };

  chemical_kinetics_fun f;
  chemical_kinetics_data_fun f_data;

  using T_t = std::tuple_element_t<2, T>;
  using Tyy = std::tuple_element_t<3, T>;
  using Typ = std::tuple_element_t<4, T>;
  using T_p = std::tuple_element_t<5, T>;

  std::vector<T_p> theta;
  Eigen::Matrix<Tyy, -1, 1> yy0;
  Eigen::Matrix<Typ, -1, 1> yp0;
  double t0;
  std::vector<T_t> ts;
  std::vector<double> x_r;
  std::vector<int> x_i;
  double rtol;
  double atol;
  int max_num_step;

  chemical_kinetics_dae_base()
      : theta{0.040, 1.0e4, 3.0e7},
        yy0(3),
        yp0(3),
        t0(0),
        ts{1, 10, 100, 1000},
        x_r{1.0e4, 3.0e7},
        rtol(1.e-5),
        atol(1.e-12),
        max_num_step(10000) {
    yy0 << 1.0, 0.0, 0.0;
    yp0 << -0.04, 0.04, 0.0;
  }

  std::vector<T_t>& times() { return ts; }
  Eigen::VectorXd init() {
    Eigen::VectorXd joined_init(yy0.size() + yp0.size());
    joined_init << yy0, yp0;
    return joined_init;
  }
  std::vector<double> param() { return theta; }
};

/**
 * Inheriting base type, various fixtures differs by the type of DAE
 * functor used in <code>apply_solver</code> calls, intended for
 * different kind of tests.
 *
 */
template <typename T>
struct chemical_kinetics_test
    : public chemical_kinetics_dae_base<T>,
      public ODETestFixture<chemical_kinetics_test<T>> {
  chemical_kinetics_test() : chemical_kinetics_dae_base<T>() {}

  auto apply_solver() {
    std::tuple_element_t<0, T> sol;
    return sol(this->f, this->yy0, this->yp0, this->t0, this->ts, nullptr,
               this->theta, this->x_r, this->x_i);
  }

  template <typename T1, typename T2>
  auto apply_solver(T1&& init, T2&& theta_in) {
    const int n = init.size() / 2;
    std::tuple_element_t<0, T> sol;
    return sol(this->f, init.head(n), init.tail(n), this->t0, this->ts, nullptr,
               theta_in, this->x_r, this->x_i);
  }

  auto apply_solver_tol() {
    std::tuple_element_t<1, T> sol;
    return sol(this->f, this->yy0, this->yp0, this->t0, this->ts, this->rtol,
               this->atol, this->max_num_step, nullptr, this->theta, this->x_r,
               this->x_i);
  }

  /**
   * Test against IDAS example "idasRoberts_FSA_dns" output.
   *
   * @param t0_in initial time
   */
  void test_value(double t0_in) {
    this->t0 = t0_in;
    this->ts.resize(3);
    this->ts[0] = this->t0 + 0.4;
    this->ts[1] = this->t0 + 4.0;
    this->ts[2] = this->t0 + 40.0;

    this->rtol = 1e-4;
    this->atol = 1e-8;
    this->max_num_step = 1000;
    {
      auto yy = apply_solver_tol();
      EXPECT_NEAR(7.1583e-01, stan::math::value_of(yy[2][0]), 2e-5);
      EXPECT_NEAR(9.1855e-06, stan::math::value_of(yy[2][1]), 1e-9);
      EXPECT_NEAR(2.8416e-01, stan::math::value_of(yy[2][2]), 2e-5);
    }

    this->rtol = 1e-8;
    this->atol = 1e-12;
    this->max_num_step = 100000;
    {
      auto yy = apply_solver_tol();
      EXPECT_NEAR(7.1583e-01, stan::math::value_of(yy[2][0]), 5e-6);
      EXPECT_NEAR(9.1855e-06, stan::math::value_of(yy[2][1]), 1e-9);
      EXPECT_NEAR(2.8416e-01, stan::math::value_of(yy[2][2]), 5e-6);
    }
  }

  /**
   * Test against IDAS example "idasRoberts_FSA_dns" output.
   *
   * @param t0_in initial time
   */
  void test_sens(double t0_in) {
    this->t0 = t0_in;
    this->ts.resize(3);
    this->ts[0] = this->t0 + 0.4;
    this->ts[1] = this->t0 + 4.0;
    this->ts[2] = this->t0 + 40.0;

    this->rtol = 1e-8;
    this->atol = 1e-12;
    this->max_num_step = 100000;
    auto yy = apply_solver_tol();
    std::vector<double> g;

    yy[2][0].grad(this->theta, g);
    EXPECT_NEAR(-4.2476e+00, g[0], 5e-5);
    EXPECT_NEAR(1.3731e-05, g[1], 1e-9);
    EXPECT_NEAR(-2.2884e-09, g[2], 1e-13);
    stan::math::set_zero_all_adjoints();

    yy[2][1].grad(this->theta, g);
    EXPECT_NEAR(4.5912e-05, g[0], 1e-8);
    EXPECT_NEAR(-2.3572e-10, g[1], 1e-12);
    EXPECT_NEAR(-1.1381e-13, g[2], 1e-15);
    stan::math::set_zero_all_adjoints();

    yy[2][2].grad(this->theta, g);
    EXPECT_NEAR(4.2475e+00, g[0], 5e-5);
    EXPECT_NEAR(-1.3731e-05, g[1], 1e-9);
    EXPECT_NEAR(2.2885e-09, g[2], 1e-13);
    stan::math::set_zero_all_adjoints();
  }

  void test_bad() {
    const auto yy0_(this->yy0);
    const auto yp0_(this->yp0);

    this->yy0.resize(0);
    EXPECT_THROW_MSG(apply_solver(), std::invalid_argument,
                     "initial state has size 0");
    this->yy0 = yy0_;

    this->yp0.resize(0);
    EXPECT_THROW_MSG(apply_solver(), std::invalid_argument,
                     "initial state derivative has size 0");
    this->yp0 = yp0_;

    const auto t0_ = this->t0;
    this->t0 = 2;
    EXPECT_THROW_MSG(apply_solver(), std::domain_error,
                     "initial time is 2, but must be less than");
    this->t0 = t0_;

    const auto ts_ = this->ts;
    this->ts.resize(0);
    EXPECT_THROW_MSG(apply_solver(), std::invalid_argument, "times has size 0");
    this->ts = ts_;

    this->ts.resize(2);
    this->ts[0] = 3.0;
    this->ts[1] = 1.0;
    EXPECT_THROW_MSG(apply_solver(), std::domain_error,
                     "times is not a valid sorted vector");
    this->ts = ts_;

    const double rtol_ = this->rtol;
    this->rtol = -1;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "relative_tolerance");
    this->rtol = rtol_;

    const double atol_ = this->atol;
    this->atol = -1;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "absolute_tolerance");
    this->atol = atol_;

    const auto theta_ = this->theta;
    const auto x_r_ = this->x_r;
    const auto x_i_ = this->x_i;

    // NaN errors
    double nan = std::numeric_limits<double>::quiet_NaN();
    std::stringstream expected_is_nan;
    expected_is_nan << "is " << nan;

    this->yy0[0] = nan;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial state");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_nan.str());
    this->yy0 = yy0_;

    this->yp0[0] = nan;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "initial state derivative");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_nan.str());
    this->yp0 = yp0_;

    this->t0 = nan;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial time");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_nan.str());
    this->t0 = t0_;

    this->ts[0] = nan;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "times");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_nan.str());
    this->ts = ts_;

    this->theta[0] = nan;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "DAE parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_nan.str());
    this->theta = theta_;

    this->x_r.push_back(nan);
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "DAE parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_nan.str());
    this->x_r = x_r_;

    // inf test
    std::stringstream expected_is_inf;
    expected_is_inf << "is " << std::numeric_limits<double>::infinity();
    std::stringstream expected_is_neg_inf;
    expected_is_neg_inf << "is " << -std::numeric_limits<double>::infinity();
    double inf = std::numeric_limits<double>::infinity();

    this->yy0[0] = inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial state");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_inf.str());
    this->yy0 = yy0_;

    this->yp0[0] = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "initial state derivative");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_neg_inf.str());
    this->yp0 = yp0_;

    this->t0 = inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial time");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_inf.str());
    this->t0 = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "initial time");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_neg_inf.str());
    this->t0 = t0_;

    this->ts.back() = inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "times");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_inf.str());
    this->ts.back() = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error, "times");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_neg_inf.str());
    this->ts = ts_;

    this->theta[0] = inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "DAE parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_inf.str());
    this->theta[0] = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "DAE parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_neg_inf.str());
    this->theta = theta_;

    this->x_r = std::vector<double>{inf};
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "DAE parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_inf.str());
    this->x_r[0] = -inf;
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     "DAE parameters and data");
    EXPECT_THROW_MSG(apply_solver_tol(), std::domain_error,
                     expected_is_neg_inf.str());
    this->x_r = x_r_;
  }
};

template <typename T>
struct chemical_kinetics_data_test
    : public chemical_kinetics_dae_base<T>,
      public ODETestFixture<chemical_kinetics_data_test<T>> {
  chemical_kinetics_data_test() : chemical_kinetics_dae_base<T>() {}

  auto apply_solver() {
    std::tuple_element_t<0, T> sol;
    return sol(this->f_data, this->yy0, this->yp0, this->t0, this->ts, nullptr,
               this->theta, this->x_r, this->x_i);
  }

  template <typename T1, typename T2>
  auto apply_solver(T1&& init, T2&& theta_in) {
    const int n = init.size() / 2;
    std::tuple_element_t<0, T> sol;
    return sol(this->f_data, init.head(n), init.tail(n), this->t0, this->ts,
               nullptr, theta_in, this->x_r, this->x_i);
  }

  auto apply_solver_tol() {
    std::tuple_element_t<1, T> sol;
    return sol(this->f_data, this->yy0, this->yp0, this->t0, this->ts,
               this->rtol, this->atol, this->max_num_step, nullptr, this->theta,
               this->x_r, this->x_i);
  }

  void test_value(double t0_in) {
    this->t0 = t0_in;
    this->ts.resize(3);
    this->ts[0] = this->t0 + 0.4;
    this->ts[1] = this->t0 + 4.0;
    this->ts[2] = this->t0 + 40.0;

    this->rtol = 1e-4;
    this->atol = 1e-8;
    this->max_num_step = 1000;
    auto yy = apply_solver_tol();
    EXPECT_NEAR(7.1583e-01, stan::math::value_of(yy[2][0]), 2e-5);
    EXPECT_NEAR(9.1855e-06, stan::math::value_of(yy[2][1]), 1e-9);
    EXPECT_NEAR(2.8416e-01, stan::math::value_of(yy[2][2]), 2e-5);
  }
};

#endif
