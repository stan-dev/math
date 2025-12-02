#include <test/unit/pretty_print_types.hpp>
#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <sstream>
#include <vector>

namespace {

class laplace_marginal_neg_binomial_log_lpmf : public LaplaceAdTest {};

TEST_P(laplace_marginal_neg_binomial_log_lpmf, phi_dim_2) {
  using stan::math::laplace_marginal_neg_binomial_2_log_lpmf;
  using stan::math::laplace_marginal_tol_neg_binomial_2_log_lpmf;
  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;

  constexpr double alpha_dbl = 1.6;
  constexpr double rho_dbl = 0.45;
  constexpr int dim_theta = 2;
  Eigen::VectorXd theta_0{{0, 0}};

  std::vector<Eigen::VectorXd> x(dim_theta);
  Eigen::VectorXd x_0{{0.05100797, 0.16086164}};
  Eigen::VectorXd x_1{{-0.59823393, 0.98701425}};
  x[0] = x_0;
  x[1] = x_1;
  std::vector<int> y{1, 0};
  std::vector<int> y_index{1, 2};
  constexpr double eta_dbl = 100;
  const auto [solver_num, hessian_block_size, max_steps_line_search]
      = GetParam();
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);
  LAPLACE_SKIP_ZERO_STEPS(max_steps_line_search);

  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 1000;
  constexpr stan::test::ad_tolerances tols{
      stan::test::ad_gradient_tols{1e-8, 1e-2}};
  auto f = [&](auto&& alpha, auto&& rho, auto&& eta) {
    try {
      return laplace_marginal_tol_neg_binomial_2_log_lpmf(
          y, y_index, eta, 0, stan::math::test::squared_kernel_functor{},
          std::forward_as_tuple(x, alpha, rho), theta_0, tolerance,
          max_num_steps, hessian_block_size, solver_num, max_steps_line_search,
          nullptr);
    } catch (const std::exception& e) {
      std::stringstream fail_msg;
      using stan::math::test::test_type_name;
      fail_msg << "Exception thrown with alpha("
               << test_type_name<decltype(alpha)>() << ")=" << alpha << ", rho("
               << test_type_name<decltype(rho)>() << ")=" << rho << ", eta("
               << test_type_name<decltype(eta)>() << ")=" << eta << ". ";
      ADD_FAILURE() << fail_msg.str() << "\n Error message: " << e.what();
      throw;
    }
  };
  stan::test::expect_ad<true>(tols, f, alpha_dbl, rho_dbl, eta_dbl);
}

LAPLACE_INSTANTIATE_TEST_SUITE_P(laplace_marginal_neg_binomial_log_lpmf);

TEST_P(laplace_disease_map_test, laplace_marginal_neg_binomial_2_log_lpmf) {
  using stan::is_var_v;
  using stan::math::laplace_marginal_neg_binomial_2_log_lpmf;
  using stan::math::laplace_marginal_tol_neg_binomial_2_log_lpmf;
  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;
  const auto [solver_num, hessian_block_size, max_steps_line_search]
      = GetParam();
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);
  LAPLACE_SKIP_ZERO_STEPS(max_steps_line_search);
  constexpr double eta = 1;

  double marginal_density = laplace_marginal_neg_binomial_2_log_lpmf(
      y, y_index, eta, mean, stan::math::test::sqr_exp_kernel_functor(),
      std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1)), nullptr);

  // TODO(charlesm93): get benchmark from GPStuff or another software.
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 1000;
  auto smoke = [&](auto&& alpha, auto&& rho, auto&& eta_arg) {
    return laplace_marginal_tol_neg_binomial_2_log_lpmf(
        y, y_index, eta_arg, mean, stan::math::test::sqr_exp_kernel_functor{},
        std::forward_as_tuple(x, alpha, rho), theta_0, tolerance, max_num_steps,
        hessian_block_size, solver_num, max_steps_line_search, nullptr);
  };
  smoke(phi_dbl[0], phi_dbl[1], eta);
  auto f = [&](auto&& alpha, auto&& rho, auto&& eta_arg) {
    try {
      return laplace_marginal_tol_neg_binomial_2_log_lpmf(
          y, y_index, eta_arg, mean, stan::math::test::sqr_exp_kernel_functor{},
          std::forward_as_tuple(x, alpha, rho), theta_0, tolerance,
          max_num_steps, hessian_block_size, solver_num, max_steps_line_search,
          nullptr);
    } catch (const std::exception& e) {
      std::stringstream fail_msg;
      using stan::math::test::test_type_name;
      fail_msg << "Exception thrown with alpha("
               << test_type_name<decltype(alpha)>() << ")=" << alpha << ", rho("
               << test_type_name<decltype(rho)>() << ")=" << rho << ", eta("
               << test_type_name<decltype(eta_arg)>() << ")=" << eta_arg
               << ". ";
      ADD_FAILURE() << fail_msg.str() << "\n Error message: " << e.what();
      throw;
    }
  };
  stan::test::expect_ad<true>(f, phi_dbl[0], phi_dbl[1], eta);
}

LAPLACE_INSTANTIATE_TEST_SUITE_P(laplace_disease_map_test);

}  // namespace
