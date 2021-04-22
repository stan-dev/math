
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
  inline auto
  operator()(const T0& t, const T1& y,
             std::ostream* msgs, const T2& ka, const T3& ke, const T4& k12,
             const T5& k21, const T6& kin, const T7& kout,
             const T8& ea50) const {
    Eigen::Matrix<stan::return_type_t<T1, T2, T3, T4, T5, T6, T7, T8>,
                  Eigen::Dynamic, 1>
        dydt(4);

    dydt << -ka * y(0), +ka * y(0) - ke * y(1) - k12 * y(1) + k21 * y(2),
        +k12 * y(1) - k21 * y(2),
        +kin - kout * (1.0 - y(1) / (y(1) + ea50)) * y(3);

    return dydt;
  }
};

void run_benchmark(int adjoint_integrator) {
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

  double log_sigma = 10.0;
  double abs_tol = 1e-8;
  double abs_tol_B = abs_tol * 100.0;
  double abs_tol_QB = abs_tol_B * 10.0;
  double rel_tol = 1e-6;
  int steps_checkpoint = 100;
  int max_num_steps = 1000000;

  for (std::size_t i = 0; i != 2; i++) {
    stan::math::nested_rev_autodiff nested;

    var CL = lognormal_rng(true_CL, log_sigma, base_rng);
    var Q = lognormal_rng(true_Q, log_sigma, base_rng);
    var V1 = lognormal_rng(true_V1, log_sigma, base_rng);
    var V2 = lognormal_rng(true_V2, log_sigma, base_rng);
    var ka = lognormal_rng(true_ka, log_sigma, base_rng);
    var pd0 = lognormal_rng(true_pd0, log_sigma, base_rng);
    var kin = lognormal_rng(true_kin, log_sigma, base_rng);
    var kout = kin / pd0;
    var ec50 = lognormal_rng(true_ec50, log_sigma, base_rng);
    var ea50 = ec50 * V1;

    var ke = CL / V1;
    var k12 = Q / V2;
    var k21 = k12 * V1 / V2;

    Eigen::Matrix<var, Eigen::Dynamic, 1> y0(4);
    y0 << 4.0, 0.0, 0.0, pd0;

    double t0 = 0.0;

    try {
      if (adjoint_integrator) {
        const int N = y0.size();
        std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> y
            = ode_adjoint_tol_ctl(
                ode, y0, t0, ts, rel_tol, Eigen::VectorXd::Constant(N, abs_tol),
                rel_tol, Eigen::VectorXd::Constant(N, abs_tol_B), rel_tol,
                abs_tol_QB, max_num_steps, steps_checkpoint, 1, 2, 2, nullptr,
                ka, ke, k12, k21, kin, kout, ea50);

        stan::math::grad();
      } else {
        std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> y
            = ode_bdf_tol(ode, y0, t0, ts, rel_tol, abs_tol_QB, max_num_steps,
                          nullptr, ka, ke, k12, k21, kin, kout, ea50);

        stan::math::grad();
      }
    } catch (std::exception& exc) {
      std::cout << "oops, keep going please!" << std::endl;
      std::cerr << exc.what() << std::endl;
    }
  }
  stan::math::recover_memory();
}

TEST(StanMathOdeBench, bdf) { run_benchmark(0); }

TEST(StanMathOdeBench, bdf_adjoint) { run_benchmark(1); }
