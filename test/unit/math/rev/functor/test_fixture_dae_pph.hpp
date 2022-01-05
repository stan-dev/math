#ifndef STAN_MATH_TEST_FIXTURE_DAE_PREY_PREDATOR_HARVEST_HPP
#define STAN_MATH_TEST_FIXTURE_DAE_PREY_PREDATOR_HARVEST_HPP

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
 * A prey-predator-harvest model. Its stability is controled by the
 * parameter mu
 */
template <typename T>
struct pph_dae_base {
  struct func {
    template <typename T0, typename Tyy, typename Typ, typename Tpar>
    inline Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> operator()(
        const T0& t, const Eigen::Matrix<Tyy, -1, 1>& yy,
        const Eigen::Matrix<Typ, -1, 1>& yp, std::ostream* msgs,
        const Tpar& mu) const {
      if (yy.size() != 3 || yp.size() != 3)
        throw std::domain_error(
            "this function was called with inconsistent state");

      Eigen::Matrix<stan::return_type_t<Tyy, Typ, Tpar>, -1, 1> res(3);

      auto yy1 = yy(0);
      auto yy2 = yy(1);
      auto yy3 = yy(2);

      auto yp1 = yp(0);
      auto yp2 = yp(1);

      const double r1 = 1.0;
      const double r2 = 3.0;
      const double p = 2.0;

      res[0] = yp1 - yy1 * (r1 - yy2);
      res[1] = yp2 - yy2 * (r2 - yy2 / yy1 - yy3);
      res[2] = yy3 * (p * yy2 - 1.0) - mu;

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

  pph_dae_base()
      : theta(0.25),
        yy0(3),
        yp0(3),
        t0(0),
        ts{0.1, 0.5, 1.0},
        rtol(1.e-5),
        atol(1.e-8),
        max_num_step(10000) {
    yy0 << 0.5, 1.0, stan::math::value_of(theta);
    yp0 << 0.0, 1.0 - stan::math::value_of(yy0[2]), 0.0;
  }

  std::vector<T_t>& times() { return ts; }
  Eigen::Matrix<stan::return_type_t<Tyy, Typ>, -1, 1> init() {
    Eigen::Matrix<stan::return_type_t<Tyy, Typ>, -1, 1> joined_init(
        yy0.size() + yp0.size());
    joined_init << stan::math::value_of(yy0), stan::math::value_of(yp0);
    return joined_init;
  }
  T_p param() { return theta; }
};

template <typename T>
struct pph_dae_test : public pph_dae_base<T>,
                      public ODETestFixture<pph_dae_test<T>> {
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
