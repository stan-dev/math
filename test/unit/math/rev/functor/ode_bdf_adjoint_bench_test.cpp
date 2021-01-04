
#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <stan/math/rev/functor/gradient.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>

using stan::math::lognormal_rng;
using stan::math::var;

struct pkpd_rhs {
  template <typename T0, typename T1, typename T2, typename T3, typename T4,
            typename T5, typename T6, typename T7, typename T8>
  inline Eigen::Matrix<stan::return_type_t<T1, T2, T3, T4, T5, T6, T7, T8>,
                       Eigen::Dynamic, 1>
  operator()(const T0& t, const Eigen::Matrix<T1, Eigen::Dynamic, 1>& y,
             std::ostream* msgs, const T2& ka, const T3& ke, const T4& k12,
             const T5& k21, const T6& kin, const T7& kout,
             const T8& ea50) const {
    Eigen::Matrix<stan::return_type_t<T1, T2, T3, T4, T5, T6, T7, T8>,
                  Eigen::Dynamic, 1>
        dydt(5);

    dydt << -ka * y(0), +ka * y(0) - ke * y(1) - k12 * y(1) + k21 * y(2),
        +k12 * y(1) - k21 * y(2),
        +kin - kout * (1.0 - y(1) / (y(1) + ea50)) * y(3),
        +kin - kout * (1.0 - y(2) / (y(2) + ea50)) * y(4);  // pseudo pd

    return dydt;
  }
};

TEST(StanMathOdeBench, bdf) {
  double true_CL = 8.0;
  double true_Q = 18.0;
  double true_V1 = 10.0;
  double true_V2 = 14.0;
  double true_ka = log(2.0) / 2.0;
  double true_pd0 = 1.0;
  double true_kin = 4.0;
  double true_kout = true_kin / true_pd0;
  double true_ec50 = 0.01;
  std::vector<double> ts{1., 2.0, 4.0, 8.0, 16.0, 32.0, 64.0};
  std::size_t ts_size = ts.size();

  pkpd_rhs ode;

  boost::ecuyer1988 base_rng(563646);

  for (std::size_t i = 0; i != 200; i++) {
    stan::math::nested_rev_autodiff nested;

    var CL = lognormal_rng(true_CL, 4.0, base_rng);
    var Q = lognormal_rng(true_Q, 4.0, base_rng);
    var V1 = lognormal_rng(true_V1, 2.0, base_rng);
    var V2 = lognormal_rng(true_V2, 2.0, base_rng);
    var ka = lognormal_rng(true_ka, 4.0, base_rng);
    var pd0 = lognormal_rng(true_pd0, 2.0, base_rng);
    var kin = lognormal_rng(true_kin, 4.0, base_rng);
    var kout = kin / pd0;
    var ec50 = lognormal_rng(true_ec50, 2.0, base_rng);
    var ea50 = ec50 * V1;

    var ke = CL / V1;
    var k12 = Q / V2;
    var k21 = k12 * V1 / V2;

    Eigen::Matrix<var, Eigen::Dynamic, 1> y0(5);
    y0 << 4.0, 0.0, 0.0, pd0, pd0 * 2.0;

    double t0 = 0.0;

    std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> y
        = ode_bdf_tol(ode, y0, t0, ts, 1E-8, 1E-8, 10000, nullptr, ka, ke, k12,
                      k21, kin, kout, ea50);

    stan::math::grad();
  }

  stan::math::recover_memory();
}

TEST(StanMathOdeBench, bdf_adjoint) {
  double true_CL = 8.0;
  double true_Q = 18.0;
  double true_V1 = 10.0;
  double true_V2 = 14.0;
  double true_ka = log(2.0) / 2.0;
  double true_pd0 = 1.0;
  double true_kin = 4.0;
  double true_kout = true_kin / true_pd0;
  double true_ec50 = 0.01;
  std::vector<double> ts{1., 2.0, 4.0, 8.0, 16.0, 32.0, 64.0};
  std::size_t ts_size = ts.size();

  pkpd_rhs ode;

  boost::ecuyer1988 base_rng(563646);

  for (std::size_t i = 0; i != 200; i++) {
    stan::math::nested_rev_autodiff nested;

    var CL = lognormal_rng(true_CL, 4.0, base_rng);
    var Q = lognormal_rng(true_Q, 4.0, base_rng);
    var V1 = lognormal_rng(true_V1, 2.0, base_rng);
    var V2 = lognormal_rng(true_V2, 2.0, base_rng);
    var ka = lognormal_rng(true_ka, 4.0, base_rng);
    var pd0 = lognormal_rng(true_pd0, 2.0, base_rng);
    var kin = lognormal_rng(true_kin, 4.0, base_rng);
    var kout = kin / pd0;
    var ec50 = lognormal_rng(true_ec50, 2.0, base_rng);
    var ea50 = ec50 * V1;

    var ke = CL / V1;
    var k12 = Q / V2;
    var k21 = k12 * V1 / V2;

    Eigen::Matrix<var, Eigen::Dynamic, 1> y0(5);
    y0 << 4.0, 0.0, 0.0, pd0, pd0 * 2.0;

    double t0 = 0.0;

    std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> y
        = ode_bdf_adjoint_tol(ode, y0, t0, ts, 1E-8, 1E-8, 10000, nullptr, ka,
                              ke, k12, k21, kin, kout, ea50);

    stan::math::grad();
  }
  stan::math::recover_memory();
}
