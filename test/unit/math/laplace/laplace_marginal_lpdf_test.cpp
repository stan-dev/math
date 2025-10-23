#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <test/unit/math/laplace/aki_synth_data/x1.hpp>
#include <test/unit/math/laplace/motorcycle_gp/x_vec.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>

struct poisson_log_likelihood2 {
  template <typename Theta>
  auto operator()(const Theta& theta, const std::vector<int>& delta_int,
                  std::ostream* pstream) const {
    return stan::math::poisson_log_lpmf(delta_int, theta);
  }
};

TEST(laplace, poisson_log_phi_dim_2) {
  using stan::math::laplace_marginal;
  using stan::math::laplace_marginal_tol;
  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;
  // logger->current_test_name_ = "poisson_log_phi_dim_2";
  constexpr int dim_phi = 2;
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_dbl{{1.6, 0.45}};

  constexpr int dim_theta = 2;
  Eigen::VectorXd theta_0(dim_theta);
  theta_0 << 0, 0;

  std::vector<Eigen::VectorXd> x(dim_theta);
  Eigen::VectorXd x_0{{0.05100797, 0.16086164}};
  Eigen::VectorXd x_1{{-0.59823393, 0.98701425}};
  x[0] = x_0;
  x[1] = x_1;

  Eigen::VectorXd y_dummy;

  std::vector<int> n_samples = {1, 1};
  std::vector<int> sums = {1, 0};

  double target = laplace_marginal<false>(
      poisson_log_likelihood2{}, std::forward_as_tuple(sums),
      stan::math::test::squared_kernel_functor{},
      std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1)), nullptr);

  // TODO(Charles): benchmark target against gpstuff.
  constexpr double tol = 1e-4;
  EXPECT_NEAR(-2.53056, value_of(target), tol);

  // Test with optional arguments
  {
    constexpr double tolerance = 1e-12;
    constexpr int max_num_steps = 100;
    constexpr int hessian_block_size = 1;
    constexpr int solver = 1;
    constexpr int max_steps_line_search = 10;

    target = laplace_marginal_tol<false>(
        poisson_log_likelihood2{}, std::forward_as_tuple(sums),
        stan::math::test::squared_kernel_functor{},
        std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1)), theta_0, tolerance,
        max_num_steps, hessian_block_size, solver, max_steps_line_search,
        nullptr);
    EXPECT_NEAR(-2.53056, value_of(target), tol);
  }

  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 100;
  using stan::is_var_v;
  using stan::scalar_type_t;
  constexpr stan::test::ad_tolerances tols{stan::test::ad_gradient_tols{1e-8, 1e-4}};
//  tols.gradient_grad_ = 1e-3;
  stan::math::test::run_solver_grid(
      [&](int solver_num, int hessian_block_size, int max_steps_line_search,
          auto&& theta_0) {
        auto f = [&](auto&& x_v, auto&& alpha, auto&& rho) {
          return laplace_marginal_tol<false>(
              poisson_log_likelihood2{}, std::forward_as_tuple(sums),
              stan::math::test::squared_kernel_functor{},
              std::forward_as_tuple(x_v, alpha, rho), theta_0, tolerance,
              max_num_steps, hessian_block_size, solver_num,
              max_steps_line_search, nullptr);
        };
        stan::test::expect_ad<true>(tols, f, x, phi_dbl[0], phi_dbl[1]);
      },
      theta_0);
}

struct poisson_log_exposure_likelihood {
  template <typename Theta, typename YEVec>
  auto operator()(const Theta& theta, YEVec&& ye,
                  const std::vector<int>& delta_int,
                  std::ostream* pstream) const {
    return stan::math::poisson_log_lpmf(
        delta_int, stan::math::add(theta, stan::math::log(ye)));
  }
};

TEST_F(laplace_disease_map_test, laplace_marginal) {

  using stan::math::laplace_marginal;
  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::laplace_marginal_tol;
  using stan::math::value_of;
  using stan::math::var;

  {
    double marginal_density = laplace_marginal<false>(
        poisson_log_exposure_likelihood{}, std::forward_as_tuple(ye, y),
        stan::math::test::sqr_exp_kernel_functor{},
        std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1)), nullptr);

    constexpr double tol = 6e-4;
    // Benchmark from GPStuff.
    EXPECT_NEAR(-2866.88, value_of(marginal_density), tol);
  }
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 100;
  stan::math::test::run_solver_grid(
      [&](int solver_num, int hessian_block_size, int max_steps_line_search,
          auto&& theta_0) {
        auto f = [&](auto&& alpha, auto&& rho) {
          return laplace_marginal_tol<false>(
              poisson_log_exposure_likelihood{}, std::forward_as_tuple(ye, y),
              stan::math::test::sqr_exp_kernel_functor{},
              std::forward_as_tuple(x, alpha, rho), theta_0, tolerance,
              max_num_steps, hessian_block_size, solver_num,
              max_steps_line_search, nullptr);
        };
        stan::test::expect_ad<true>(f, phi_dbl[0], phi_dbl[1]);
      },
      theta_0);
}

struct bernoulli_logit_likelihood {
  template <typename Theta>
  auto operator()(const Theta& theta, const std::vector<int>& delta_int,
                  std::ostream* pstream) const {
    return stan::math::bernoulli_logit_lpmf(delta_int, theta);
  }
};

TEST(laplace, bernoulli_logit_phi_dim500) {

  using stan::math::laplace_marginal;
  using stan::math::laplace_marginal_tol;
  using stan::math::to_vector;
  // logger->current_test_name_ = "bernoulli_logit_phi_dim500";
  constexpr int dim_theta = 500;
  constexpr int n_observations = 500;
  auto x1 = stan::test::laplace::x1;
  auto x2 = stan::test::laplace::x2;
  auto y = stan::test::laplace::y;

  constexpr int dim_x = 2;
  std::vector<Eigen::VectorXd> x(dim_theta);
  for (int i = 0; i < dim_theta; i++) {
    Eigen::VectorXd coordinate(dim_x);
    coordinate << x1[i], x2[i];
    x[i] = coordinate;
  }
  Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);
  Eigen::VectorXd delta_L;
  std::vector<double> delta;
  constexpr int dim_phi = 2;
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_dbl{{1.6, 1}};

  double target = laplace_marginal<false>(
      bernoulli_logit_likelihood{}, std::forward_as_tuple(y),
      stan::math::test::sqr_exp_kernel_functor{},
      std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1)), nullptr);

  constexpr double tol = 3e-4;
  // Benchmark against gpstuff.
  EXPECT_NEAR(-195.368, target, tol);
  // All fail for ad check with relative tolerance ~0.002
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 100;
  constexpr stan::test::ad_tolerances tols;
  stan::math::test::run_solver_grid(
      [&](int solver_num, int hessian_block_size, int max_steps_line_search,
          auto&& theta_0) {
        auto f = [&](auto&& alpha, auto&& rho) {
          return laplace_marginal_tol<false>(
              bernoulli_logit_likelihood{}, std::forward_as_tuple(y),
              stan::math::test::sqr_exp_kernel_functor{},
              std::forward_as_tuple(x, alpha, rho), theta_0, tolerance,
              max_num_steps, hessian_block_size, solver_num,
              max_steps_line_search, nullptr);
        };
        stan::test::expect_ad<true>(tols, f, phi_dbl[0], phi_dbl[1]);
      },
      theta_0);
}

