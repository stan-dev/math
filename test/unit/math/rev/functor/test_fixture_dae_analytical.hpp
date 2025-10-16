#ifndef STAN_MATH_TEST_FIXTURE_DAE_ANALYTICAL_HPP
#define STAN_MATH_TEST_FIXTURE_DAE_ANALYTICAL_HPP

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
struct analytical_dae_base {
  struct func {
    template <typename T0, typename Tyy, typename Typ, typename Tpar>
    inline Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> operator()(
        const T0& t_in, const Eigen::Matrix<Tyy, -1, 1>& yy,
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

      res[0] = yp1 - theta * yy3;
      res[1] = yy2 * (yy2 - 1.0);
      res[2] = t_in - yy1 * yy2 - yy3 * (1.0 - yy2);

      return res;
    }
  };

  func f;

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

  analytical_dae_base()
      : theta(2.0),
        yy0(3),
        yp0(3),
        t0(0),
        ts{1, 10, 100},
        rtol(1.e-5),
        atol(1.e-12),
        max_num_step(10000) {
    yy0 << 1.0, 0.0, 0.0;
    yp0 << 0.0, 0.77, 0.77;
  }

  std::vector<T_t>& times() { return ts; }
  Eigen::VectorXd init() {
    Eigen::VectorXd joined_init(yy0.size() + yp0.size());
    joined_init << yy0, yp0;
    return joined_init;
  }
  std::vector<double> param() { return theta; }
};

template <typename T>
struct analytical_dae_dv_functor : public analytical_dae_base<T> {
  analytical_dae_dv_functor() : analytical_dae_base<T>() {}

  template <typename Tx>
  Eigen::Matrix<Tx, -1, 1> operator()(Eigen::Matrix<Tx, -1, 1>& x) const {
    std::tuple_element_t<0, T> sol;
    auto ys
        = sol(this->f, this->yy0, this->yp0, this->t0, this->ts, nullptr, x(0));
    return ys[0];
  }
};

template <typename T>
struct analytical_dae_vd_functor : public analytical_dae_base<T> {
  analytical_dae_vd_functor() : analytical_dae_base<T>() {}

  template <typename Tx>
  Eigen::Matrix<Tx, -1, 1> operator()(Eigen::Matrix<Tx, -1, 1>& x) const {
    std::tuple_element_t<0, T> sol;
    Eigen::Matrix<Tx, -1, 1> yy0_tx = x;
    auto ys = sol(this->f, yy0_tx, this->yp0, this->t0, this->ts, nullptr,
                  this->theta);
    return ys[0];
  }
};

template <typename T>
struct analytical_dae_test : public analytical_dae_base<T>,
                             public ODETestFixture<analytical_dae_test<T>> {
  analytical_dae_dv_functor<T> dae_sol_dv;
  analytical_dae_vd_functor<T> dae_sol_vd;
  // analytical_dae_vv_functor<T> dae_sol_vv;

  auto analy_sol_functor() {
    auto f = [&](double t, double k) {
      Eigen::VectorXd y(3);
      y << 0.5 * k * t * t + stan::math::value_of(this->yy0[0]), 0, t;
      return y;
    };
    return f;
  }

  auto analy_grad_theta_functor() {
    auto f = [&](double t, double k) {
      Eigen::MatrixXd dy(1, 3);
      dy(0, 0) = 0.5 * t * t;
      dy(0, 1) = 0;
      dy(0, 2) = 0;
      return dy;
    };
    return f;
  }

  auto analy_grad_yy0_functor() {
    auto f = [&](double t, double k) {
      Eigen::MatrixXd dy(3, 3);

      // dyy[0]/d(yy0)
      dy(0, 0) = 1;
      dy(1, 0) = 0;
      dy(2, 0) = 0;

      // dyy[1]/d(yy0)
      dy(0, 1) = 0;
      dy(1, 1) = 0;
      dy(2, 1) = 0;

      // dyy[2]/d(yy0)
      dy(0, 2) = 0;
      dy(1, 2) = 0;
      dy(2, 2) = 0;
      return dy;
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
               theta_in);
  }

  auto apply_solver_tol() {
    std::tuple_element_t<1, T> sol;
    return sol(this->f, this->yy0, this->yp0, this->t0, this->ts, this->rtol,
               this->atol, this->max_num_step, nullptr, this->theta);
  }
};

#endif
