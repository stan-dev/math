#ifndef STAN_MATH_TEST_FIXTURE_DAE_INDEX_3_HPP
#define STAN_MATH_TEST_FIXTURE_DAE_INDEX_3_HPP

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
 * A DAE of index 3
 */
template <typename T>
struct index_3_dae_base {
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

      res[0] = yp1 + yy2 - sin(t);
      res[1] = yp2 + yy3 - sin(t);
      res[2] = yy1 - cos(t);

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

  index_3_dae_base()
      : theta(0.25),
        yy0(3),
        yp0(3),
        t0(0),
        ts{0.1, 1.0, 10.0},
        rtol(1.e-5),
        atol(1.e-8),
        max_num_step(10000000) {
    yy0 << 1.0, 0.0, -2.0;
    yp0 << 0.0, 2.0, 1.0;
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
struct index_3_dae_test : public index_3_dae_base<T>,
                          public ODETestFixture<index_3_dae_test<T>> {
  auto apply_solver() {
    std::tuple_element_t<0, T> sol;
    return sol(this->f, this->yy0, this->yp0, this->t0, this->ts, nullptr,
               this->theta);
  }

  auto apply_solver_tol() {
    std::tuple_element_t<1, T> sol;
    return sol(this->f, this->yy0, this->yp0, this->t0, this->ts, this->rtol,
               this->atol, this->max_num_step, nullptr, this->theta);
  }
};

#endif
