

#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <stan/math/rev/functor/gradient.hpp>
#include <iostream>
#include <sstream>
#include <vector>
#include <stdexcept>

using stan::math::lognormal_rng;
using stan::math::var;

struct scaling_rhs {
  template <typename T0, typename T1, typename T2, typename T3, typename T4,
            typename T5>
  inline auto
  operator()(const T0& t, const T1& y,
             std::ostream* msgs, const T2& kt,
             const T3& e50, const T4& k12,
             const T5& k21) const {
    std::size_t num_main_states = kt.size();
    std::size_t num_states = 2 * num_main_states;

    using return_t = stan::return_type_t<T1, T2, T3, T4, T5>;

    Eigen::Matrix<return_t, Eigen::Dynamic, 1> dydt(num_states);
    std::vector<return_t> ksat(num_main_states);

    for (std::size_t i = 0; i != num_main_states; ++i) {
      std::size_t m = 2 * i;  // main state
      std::size_t a = m + 1;  // auxilary state
      ksat[i] = kt[i] * y(m) / (y(m) + e50[i]);

      dydt(m) = -ksat[i] * y(m) - k12[i] * y(m) + k21[i] * y(a);
      dydt(a) = +k12[i] * y(m) - k21[i] * y(a);

      if (i != 0) {
        dydt(m) += ksat[i - 1] * y(2 * (i - 1));
      }
    }
    return dydt;
  }
};

void run_benchmark(std::size_t system_size, int adjoint_integrator) {
  scaling_rhs ode;

  std::vector<double> ts{1., 2.0, 4.0, 8.0, 16.0, 32.0, 64.0};
  std::size_t ts_size = ts.size();

  boost::ecuyer1988 base_rng(563646);

  double log_sigma = 10.0;
  double abs_tol = 1e-8;
  double abs_tol_B = abs_tol * 10.0;
  double abs_tol_QB = abs_tol_B * 10.0;
  double rel_tol = 1e-6;
  int steps_checkpoint = 100;
  int max_num_steps = 1000000;

  for (std::size_t i = 0; i != 2; i++) {
    stan::math::nested_rev_autodiff nested;

    std::vector<var> kt(system_size);
    std::vector<var> e50(system_size);
    std::vector<var> k12(system_size);
    std::vector<var> k21(system_size);

    for (std::size_t j = 0; j != system_size; ++j) {
      kt[j] = lognormal_rng(0.0, log_sigma, base_rng);
      e50[j] = lognormal_rng(0.0, log_sigma, base_rng);
      k12[j] = lognormal_rng(0.0, log_sigma, base_rng);
      k21[j] = lognormal_rng(0.0, log_sigma, base_rng);
    }

    Eigen::Matrix<var, Eigen::Dynamic, 1> y0(2 * system_size);
    for (std::size_t j = 0; j != 2 * system_size; ++j)
      y0(j) = lognormal_rng(0.0, log_sigma, base_rng);

    double t0 = 0.0;

    try {
      std::vector<Eigen::Matrix<var, Eigen::Dynamic, 1>> y;
      if (adjoint_integrator) {
        const int N = y0.size();
        y = ode_adjoint_tol_ctl(
            ode, y0, t0, ts, rel_tol, Eigen::VectorXd::Constant(N, abs_tol),
            rel_tol, Eigen::VectorXd::Constant(N, abs_tol_B), rel_tol,
            abs_tol_QB, max_num_steps, steps_checkpoint, 1, 2, 2, nullptr, kt,
            e50, k12, k21);
      } else {
        y = ode_bdf_tol(ode, y0, t0, ts, rel_tol, abs_tol_QB, max_num_steps,
                        nullptr, kt, e50, k12, k21);
      }

      // Essentially sets the adjoint for all states to 1.
      var target = stan::math::sum(y[0]);
      for (int k = 1; k < ts_size; k++)
        target += stan::math::sum(y[k]);

      target.grad();
    } catch (std::exception& exc) {
      std::cout << "oops, keep going please!" << std::endl;
      std::cerr << exc.what() << std::endl;
    }
  }
  stan::math::recover_memory();
}

TEST(StanMathOdeBench, bdf_2) { run_benchmark(2, 0); }

TEST(StanMathOdeBench, bdf_adjoint_2) { run_benchmark(2, 1); }

TEST(StanMathOdeBench, bdf_4) { run_benchmark(4, 0); }

TEST(StanMathOdeBench, bdf_adjoint_4) { run_benchmark(4, 1); }

TEST(StanMathOdeBench, bdf_8) { run_benchmark(8, 0); }

TEST(StanMathOdeBench, bdf_adjoint_8) { run_benchmark(8, 1); }

TEST(StanMathOdeBench, bdf_16) { run_benchmark(16, 0); }

TEST(StanMathOdeBench, bdf_adjoint_16) { run_benchmark(16, 1); }

TEST(StanMathOdeBench, bdf_32) { run_benchmark(32, 0); }

TEST(StanMathOdeBench, bdf_adjoint_32) { run_benchmark(32, 1); }
