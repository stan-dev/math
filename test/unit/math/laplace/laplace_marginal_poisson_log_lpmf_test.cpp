#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <vector>

TEST(laplace_marginal_poisson_log_lpmf, phi_dim_2) {
  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::laplace_marginal_tol_poisson_log_lpmf;

  using stan::math::log;
  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;

  double alpha_dbl = 1.6;
  double rho_dbl = 0.45;
  int dim_theta = 2;
  Eigen::VectorXd theta_0(dim_theta);
  theta_0 << 0, 0;

  std::vector<Eigen::VectorXd> x(dim_theta);
  Eigen::VectorXd x_0(2);
  x_0 << 0.05100797, 0.16086164;
  Eigen::VectorXd x_1(2);
  x_1 << -0.59823393, 0.98701425;
  x[0] = x_0;
  x[1] = x_1;

  std::vector<double> delta;
  std::vector<int> delta_int;

  std::vector<int> y = {1, 0};
  std::vector<int> y_index = {1, 2};

  stan::math::test::squared_kernel_functor sq_kernel;
  constexpr double tolerance = 1e-6;
  constexpr int max_num_steps = 100;

  stan::test::ad_tolerances tols;
  // tols.gradient_val_ = 1e-3;
  tols.gradient_grad_ = 1e-3;

  for (int max_steps_line_search = 0; max_steps_line_search < 4;
       ++max_steps_line_search) {
    for (int hessian_block_size = 1; hessian_block_size < 4;
         hessian_block_size++) {
      for (int solver_num = 1; solver_num < 4; solver_num++) {
        auto f = [&](auto&& alpha, auto&& rho) {
          return laplace_marginal_tol_poisson_log_lpmf(
              y, y_index, 0, sq_kernel, std::forward_as_tuple(x, alpha, rho),
              theta_0, tolerance, max_num_steps, hessian_block_size, solver_num,
              max_steps_line_search, nullptr);
        };
        stan::test::expect_ad<true>(tols, f, alpha_dbl, rho_dbl);
      }
    }
  }

  Eigen::VectorXd ye(2);
  ye << 1, 1;
  for (int max_steps_line_search = 0; max_steps_line_search < 4;
       ++max_steps_line_search) {
    for (int hessian_block_size = 1; hessian_block_size < 4;
         hessian_block_size++) {
      for (int solver_num = 1; solver_num < 4; solver_num++) {
        auto f = [&](auto&& alpha, auto&& rho) {
          return laplace_marginal_tol_poisson_log_lpmf(
              y, y_index, log(ye), sq_kernel,
              std::forward_as_tuple(x, alpha, rho), theta_0, tolerance,
              max_num_steps, hessian_block_size, solver_num,
              max_steps_line_search, nullptr);
        };
        stan::test::expect_ad<true>(tols, f, alpha_dbl, rho_dbl);
      }
    }
  }
}

TEST_F(laplace_disease_map_test, laplace_marginal_poisson_log_lpmf) {
  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::laplace_marginal_tol_poisson_log_lpmf;
  using stan::math::log;
  using stan::math::value_of;
  using stan::math::var;

  double marginal_density = laplace_marginal_poisson_log_lpmf(
      y, y_index, log(ye), stan::math::test::sqr_exp_kernel_functor(),
      std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1)), nullptr);

  double tol = 6e-4;
  // Benchmark from GPStuff.
  EXPECT_NEAR(-2866.88, marginal_density, tol);

  constexpr double tolerance = 1e-6;
  constexpr int max_num_steps = 100;
  for (int max_steps_line_search = 0; max_steps_line_search < 4;
       ++max_steps_line_search) {
    for (int hessian_block_size = 1; hessian_block_size < 4;
         hessian_block_size++) {
      for (int solver_num = 1; solver_num < 4; solver_num++) {
        auto f = [&](auto&& alpha, auto&& rho) {
          return laplace_marginal_tol_poisson_log_lpmf(
              y, y_index, log(ye), stan::math::test::sqr_exp_kernel_functor(),
              std::forward_as_tuple(x, alpha, rho), theta_0, tolerance,
              max_num_steps, hessian_block_size, solver_num,
              max_steps_line_search, nullptr);
        };
        stan::test::expect_ad<true>(f, phi_dbl[0], phi_dbl[1]);
      }
    }
  }
}

struct diag_covariance {
  template <typename T0__>
  Eigen::Matrix<stan::return_type_t<T0__>, -1, -1> operator()(
      const T0__& sigma, const int& N, std::ostream* pstream__) const {
    return stan::math::diag_matrix(
        stan::math::rep_vector(stan::math::pow(sigma, 2), N));
  }
};

TEST(laplace_marginal_poisson_log_lpmf, mean_argument) {
  // working example from
  // https://discourse.mc-stan.org/t/embedded-laplace-numerical-problem/39700
  using stan::math::laplace_marginal_poisson_log_lpmf;

  const int N = 1;
  const std::vector<int> y{153};
  const std::vector<int> y_index{1};

  Eigen::VectorXd mu(1);
  mu << 4.3;

  const double sigmaz = 2.0;

  double marginal_density = laplace_marginal_poisson_log_lpmf(
      y, y_index, mu, diag_covariance(), std::tuple<double, int>(sigmaz, N),
      nullptr);

  EXPECT_FLOAT_EQ(-6.7098737, marginal_density);
}
