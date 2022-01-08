#ifndef STAN_MATH_TEST_FIXTURE_DAE_LINEAR_HPP
#define STAN_MATH_TEST_FIXTURE_DAE_LINEAR_HPP

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

/**
 * Linear constant coef DAE in the form
 * x''(t) = z(t)
 * x(t) + z(t) = 3 sin(t)
 *
 * with IC x(pi) = 1.0, x'(pi) = 0.0
 *
 */
template <typename T>
struct linear_dae_base {
  struct func {
    template <typename T0, typename Tyy, typename Typ, typename Tpar>
    inline Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> operator()(
        const T0& t, const Eigen::Matrix<Tyy, -1, 1>& yy,
        const Eigen::Matrix<Typ, -1, 1>& yp, std::ostream* msgs,
        const Tpar& theta) const {
      if (yy.size() != 3 || yp.size() != 3)
        throw std::domain_error(
            "this function was called with inconsistent state");

      Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> res(3);

      auto yy1 = yy(0);
      auto yy2 = yy(1);
      auto yy3 = yy(2);

      auto yp1 = yp(0);
      auto yp2 = yp(1);

      res[0] = yp1 - yy2;
      res[1] = yp2 - yy3;
      res[2] = yy1 + yy3 - theta * sin(t);

      return res;
    }
  };

  struct degen_func {
    template <typename T0, typename Tyy, typename Typ, typename Tpar>
    inline Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> operator()(
        const T0& t, const Eigen::Matrix<Tyy, -1, 1>& yy,
        const Eigen::Matrix<Typ, -1, 1>& yp, std::ostream* msgs,
        const Tpar& theta) const {
      Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> res(2);
      res[0] = yp(0) - yy(1);
      res[1] = yp(1) + yy(0) - theta * sin(t);
      return res;
    }
  };

  struct ode_func {
    template <typename T0, typename Ty, typename Tpar>
    inline Eigen::Matrix<stan::return_type_t<Ty, Tpar>, -1, 1> operator()(
        const T0& t, const Eigen::Matrix<Ty, -1, 1>& y,
        const std::ostream* msgs, const Tpar& theta) const {
      Eigen::Matrix<stan::return_type_t<Ty, Tpar>, -1, 1> dydt(2);
      dydt << y[1], theta * sin(t) - y[0];
      return dydt;
    }
  };

  func f;
  degen_func f_degen;
  ode_func f_ode;

  using T_t = std::tuple_element_t<2, T>;
  using Tyy = std::tuple_element_t<3, T>;
  using Typ = std::tuple_element_t<4, T>;
  using T_p = std::tuple_element_t<5, T>;

  T_p theta;
  Eigen::Matrix<Tyy, -1, 1> yy0;
  Eigen::Matrix<Typ, -1, 1> yp0;
  double t0;
  std::vector<T_t> ts;
  double rtol;
  double atol;
  int max_num_step;

  linear_dae_base()
      : theta(3.0),
        yy0(3),
        yp0(3),
        t0(3.14159265358979323846),
        ts{2.4 * t0},
        rtol(1.e-5),
        atol(1.e-8),
        max_num_step(1000) {
    yy0 << 1.0, 0.0, -1.0;
    yp0 << 0.0, -1.0, 0.0;
  }

  std::vector<T_t>& times() { return ts; }
  Eigen::VectorXd init() {
    Eigen::VectorXd joined_init(yy0.size() + yp0.size());
    joined_init << stan::math::value_of(yy0), stan::math::value_of(yp0);
    return joined_init;
  }
  std::vector<double> param() { return {theta}; }
};

template <typename T>
struct linear_dae_functor : public linear_dae_base<T> {
  template <typename Tx>
  Eigen::Matrix<Tx, -1, 1> operator()(const Eigen::Matrix<Tx, -1, 1>& x) const {
    std::tuple_element_t<0, T> sol;
    auto ys
        = sol(this->f, this->yy0, this->yp0, this->t0, this->ts, nullptr, x[0]);
    return ys[0];
  }
};

template <typename T>
struct linear_dae_test : public linear_dae_base<T>,
                         public ODETestFixture<linear_dae_test<T>> {
  linear_dae_functor<T> dae_sol;

  /**
   * analytical solution based on theta = 3.0
   *
   */
  auto analy_sol_functor() {
    auto f = [&](double t) {
      const double pi = 3.14159265358979323846;
      Eigen::VectorXd y(3);
      y << -cos(t) + 1.5 * pi * cos(t) - 1.5 * t * cos(t) + 1.5 * sin(t),
          sin(t) - 1.5 * pi * sin(t) - 1.5 * cos(t) + 1.5 * t * sin(t)
              + 1.5 * cos(t),
          cos(t) - 1.5 * pi * cos(t) + 1.5 * t * cos(t) + 1.5 * sin(t);
      return y;
    };
    return f;
  }

  auto apply_solver() {
    std::tuple_element_t<0, T> sol;
    return sol(this->f, this->yy0, this->yp0, this->t0, this->ts, nullptr,
               this->theta);
  }

  template <typename T1, typename T2>
  auto apply_solver(T1&& init, T2&& theta_in) {
    const int n = init.size() / 2;
    std::tuple_element_t<0, T> sol;
    return sol(this->f, init.head(n), init.tail(n), this->t0, this->ts, nullptr,
               theta_in[0]);
  }

  auto apply_solver_tol() {
    std::tuple_element_t<1, T> sol;
    return sol(this->f, this->yy0, this->yp0, this->t0, this->ts, this->rtol,
               this->atol, this->max_num_step, nullptr, this->theta);
  }
};

template <typename T>
struct degenerated_dae_test : public linear_dae_base<T>,
                              public ODETestFixture<linear_dae_test<T>> {
  void test_ode_sens_y0() {
    this->ts = {4.0, 6.0};

    std::tuple_element_t<0, T> sol;

    stan::math::nested_rev_autodiff nested;
    Eigen::Matrix<stan::math::var, -1, 1> dae_y0_var(2);
    dae_y0_var << 1.0, 0.0;
    Eigen::Matrix<stan::math::var, -1, 1> ode_y0_var(2);
    ode_y0_var << 1.0, 0.0;

    auto res1 = sol(this->f_degen, dae_y0_var, this->yp0, this->t0, this->ts,
                    nullptr, this->theta);
    auto res2 = stan::math::ode_bdf(this->f_ode, ode_y0_var, this->t0, this->ts,
                                    nullptr, this->theta);

    std::vector<double> g1(2), g2(2);
    for (auto i = 0; i < 2; ++i) {
      nested.set_zero_all_adjoints();
      res1[1][i].grad();
      g1[0] = dae_y0_var[0].adj();
      g1[1] = dae_y0_var[1].adj();
      nested.set_zero_all_adjoints();
      res2[1][i].grad();
      g2[0] = ode_y0_var[0].adj();
      g2[1] = ode_y0_var[1].adj();
      EXPECT_NEAR(g1[0], g2[0], 4.e-9);
      EXPECT_NEAR(g1[1], g2[1], 4.e-9);
    }
  }

  void test_ode_sens_theta() {
    this->ts = {4.0, 6.0};
    this->yy0.resize(2);
    this->yy0[0] = 1.0;
    this->yy0[1] = 0.0;
    this->yp0.resize(2);
    this->yp0[0] = 0.0;
    this->yp0[1] = -1.0;

    std::tuple_element_t<0, T> sol;

    stan::math::nested_rev_autodiff nested;
    stan::math::var dae_theta(3.0), ode_theta(3.0);
    auto res1 = sol(this->f_degen, this->yy0, this->yp0, this->t0, this->ts,
                    nullptr, dae_theta);
    auto res2 = stan::math::ode_bdf(this->f_ode, this->yy0, this->t0, this->ts,
                                    nullptr, ode_theta);

    double g1, g2;
    for (auto i = 0; i < 2; ++i) {
      nested.set_zero_all_adjoints();
      res1[1][i].grad();
      g1 = dae_theta.adj();
      nested.set_zero_all_adjoints();
      res2[1][i].grad();
      g2 = ode_theta.adj();
      EXPECT_NEAR(g1, g2, 5.e-9);
    }
  }
};

#endif
