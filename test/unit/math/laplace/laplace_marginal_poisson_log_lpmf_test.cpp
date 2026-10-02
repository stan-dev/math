#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/pretty_print_types.hpp>

#include <gtest/gtest.h>
#include <sstream>
#include <vector>

namespace {

class laplace_marginal_poisson_log_lpmf : public LaplaceAdTest {};

TEST_P(laplace_marginal_poisson_log_lpmf, phi_dim_2) {
  constexpr int dim_theta = 2;
  const auto test_params = GetParam();
  const auto solver_num = std::get<0>(test_params);
  const auto hessian_block_size = std::get<1>(test_params);
  const auto max_steps_line_search = std::get<2>(test_params);
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);

  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::laplace_marginal_tol_poisson_log_lpmf;

  using stan::math::log;
  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;

  double alpha_dbl = 1.6;
  double rho_dbl = 0.45;
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
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 500;
  constexpr bool allow_fallthrough = true;

  constexpr stan::test::ad_tolerances tols;
  // tols.gradient_val_ = 1e-3;
  //  tols.gradient_grad_ = 1e-3;
  auto f = [&](auto&& alpha, auto&& rho) {
    try {
      return laplace_marginal_tol_poisson_log_lpmf(
          y, y_index, 0, hessian_block_size, sq_kernel,
          std::forward_as_tuple(x, alpha, rho),
          std::make_tuple(theta_0, tolerance, max_num_steps, solver_num,
                          max_steps_line_search, allow_fallthrough),
          &output_stream);
    } catch (const std::exception& e) {
      std::stringstream fail_msg;
      using stan::math::test::test_type_name;
      fail_msg << "Exception thrown with alpha("
               << test_type_name<decltype(alpha)>() << ")=" << alpha << ", rho("
               << test_type_name<decltype(rho)>() << ")=" << rho << ". ";
      ADD_FAILURE() << fail_msg.str() << "\n Error message: " << e.what();
      throw;
    }
  };
  stan::test::expect_ad<true>(tols, f, alpha_dbl, rho_dbl);
}

TEST_P(laplace_marginal_poisson_log_lpmf, log_phi_dim_2) {
  constexpr int dim_theta = 2;
  const auto test_params = GetParam();
  const auto solver_num = std::get<0>(test_params);
  const auto hessian_block_size = std::get<1>(test_params);
  const auto max_steps_line_search = std::get<2>(test_params);
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);
  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::laplace_marginal_tol_poisson_log_lpmf;

  using stan::math::log;
  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;

  double alpha_dbl = 1.6;
  double rho_dbl = 0.45;
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
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 500;

  //  stan::test::ad_tolerances tols;
  // tols.gradient_val_ = 1e-3;
  constexpr stan::test::ad_tolerances tols{
      stan::test::ad_gradient_tols{1e-8, 1e-3}};
  //  tols.gradient_grad_ = 1e-3;
  Eigen::VectorXd ye(2);
  ye << 1, 1;
  auto f = [&](auto&& alpha, auto&& rho) {
    try {
      return laplace_marginal_tol_poisson_log_lpmf(
          y, y_index, log(ye), hessian_block_size, sq_kernel,
          std::forward_as_tuple(x, alpha, rho),
          std::make_tuple(theta_0, tolerance, max_num_steps, solver_num,
                          max_steps_line_search, true),
          &output_stream);
    } catch (const std::exception& e) {
      std::stringstream fail_msg;
      using stan::math::test::test_type_name;
      fail_msg << "Exception thrown with alpha("
               << test_type_name<decltype(alpha)>() << ")=" << alpha << ", rho("
               << test_type_name<decltype(rho)>() << ")=" << rho << ". ";
      ADD_FAILURE() << fail_msg.str() << "\n Error message: " << e.what();
      throw;
    }
  };
  stan::test::expect_ad<true>(tols, f, alpha_dbl, rho_dbl);
}

struct diag_covariance {
  template <typename T0__>
  Eigen::Matrix<stan::return_type_t<T0__>, -1, -1> operator()(
      const T0__& sigma, const int& N, std::ostream* pstream__) const {
    return stan::math::diag_matrix(
        stan::math::rep_vector(stan::math::pow(sigma, 2), N));
  }
};

TEST_P(laplace_marginal_poisson_log_lpmf, mean_argument) {
  constexpr int dim_theta = 1;
  const auto test_params = GetParam();
  const auto solver_num = std::get<0>(test_params);
  const auto hessian_block_size = std::get<1>(test_params);
  const auto max_steps_line_search = std::get<2>(test_params);
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);

  // working example from
  // https://discourse.mc-stan.org/t/embedded-laplace-numerical-problem/39700
  using stan::math::laplace_marginal_tol_poisson_log_lpmf;

  const std::vector<int> y{153};
  const std::vector<int> y_index{1};
  Eigen::VectorXd theta_0(dim_theta);
  theta_0 << 0;

  Eigen::VectorXd mu(1);
  mu << 4.3;

  constexpr double sigmaz = 2.0;
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 500;
  double marginal_density = laplace_marginal_tol_poisson_log_lpmf(
      y, y_index, mu, hessian_block_size, diag_covariance{},
      std::tuple<double, int>(sigmaz, dim_theta),
      std::make_tuple(theta_0, tolerance, max_num_steps, solver_num,
                      max_steps_line_search, true),
      &output_stream);

  EXPECT_FLOAT_EQ(-6.7098737, marginal_density);
}
LAPLACE_INSTANTIATE_TEST_SUITE_P(laplace_marginal_poisson_log_lpmf);

TEST_P(laplace_disease_map_test, laplace_marginal_poisson_log_lpmf) {
  const auto test_params = GetParam();
  const auto solver_num = std::get<0>(test_params);
  const auto hessian_block_size = std::get<1>(test_params);
  const auto max_steps_line_search = std::get<2>(test_params);
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);

  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::laplace_marginal_tol_poisson_log_lpmf;
  using stan::math::log;
  using stan::math::value_of;
  using stan::math::var;
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 500;

  double marginal_density = laplace_marginal_tol_poisson_log_lpmf(
      y, y_index, log(ye), hessian_block_size,
      stan::math::test::sqr_exp_kernel_functor(),
      std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1)),
      std::make_tuple(theta_0, tolerance, max_num_steps, solver_num,
                      max_steps_line_search, true),
      &output_stream);

  double tol = 6e-4;
  // Benchmark from GPStuff.
  EXPECT_NEAR(-2866.88, marginal_density, tol);
  auto f = [&](auto&& alpha, auto&& rho) {
    try {
      return laplace_marginal_tol_poisson_log_lpmf(
          y, y_index, log(ye), hessian_block_size,
          stan::math::test::sqr_exp_kernel_functor(),
          std::forward_as_tuple(x, alpha, rho),
          std::make_tuple(theta_0, tolerance, max_num_steps, solver_num,
                          max_steps_line_search, true),
          &output_stream);
    } catch (const std::exception& e) {
      std::stringstream fail_msg;
      using stan::math::test::test_type_name;
      fail_msg << "Exception thrown with alpha("
               << test_type_name<decltype(alpha)>() << ")=" << alpha << ", rho("
               << test_type_name<decltype(rho)>() << ")=" << rho << ". ";
      ADD_FAILURE() << fail_msg.str() << "\n Error message: " << e.what();
      throw;
    }
  };
  stan::test::expect_ad<true>(f, phi_dbl[0], phi_dbl[1]);
}
LAPLACE_INSTANTIATE_TEST_SUITE_P(laplace_disease_map_test);

// Reference likelihood computed per observation, without the grouped
// sufficient-statistics shortcut used by poisson_log_likelihood.
struct poisson_log_obs_likelihood {
  template <typename Theta, typename Mean>
  auto operator()(const Theta& theta, const std::vector<int>& y,
                  const std::vector<int>& y_index, const Mean& mean,
                  std::ostream* /*pstream*/) const {
    Eigen::Matrix<stan::return_type_t<Theta, Mean>, Eigen::Dynamic, 1>
        theta_obs(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
      theta_obs(i) = theta(y_index[i] - 1) + mean(y_index[i] - 1);
    }
    return stan::math::poisson_log_lpmf(y, theta_obs);
  }
};

// More observations than latent variables: groups with several observations
// must aggregate all of them, not just the first theta.size() entries.
TEST(laplace_marginal_poisson_log_lpmf_grouped, multiple_obs_per_group) {
  using stan::math::laplace_marginal_tol;
  using stan::math::laplace_marginal_tol_poisson_log_lpmf;

  constexpr int dim_theta = 3;
  const std::vector<int> y{1, 0, 3, 2, 4, 1};
  const std::vector<int> y_index{1, 2, 3, 1, 2, 1};
  const Eigen::VectorXd mean{{0.3, -0.2, 0.1}};
  const std::vector<Eigen::VectorXd> x{
      Eigen::VectorXd{{0.05100797, 0.16086164}},
      Eigen::VectorXd{{-0.59823393, 0.98701425}},
      Eigen::VectorXd{{0.31296868, -0.68926772}}};
  const Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);
  constexpr double alpha = 1.6;
  constexpr double rho = 0.45;
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 1000;
  constexpr int hessian_block_size = 1;
  constexpr int solver = 1;
  constexpr int max_steps_line_search = 0;

  const double marginal = laplace_marginal_tol_poisson_log_lpmf(
      y, y_index, mean, hessian_block_size,
      stan::math::test::squared_kernel_functor{},
      std::forward_as_tuple(x, alpha, rho),
      std::make_tuple(theta_0, tolerance, max_num_steps, solver,
                      max_steps_line_search, true),
      nullptr);
  const double reference = laplace_marginal_tol<false>(
      poisson_log_obs_likelihood{}, std::forward_as_tuple(y, y_index, mean),
      hessian_block_size, stan::math::test::squared_kernel_functor{},
      std::forward_as_tuple(x, alpha, rho),
      std::make_tuple(theta_0, tolerance, max_num_steps, solver,
                      max_steps_line_search, true),
      nullptr);
  EXPECT_NEAR(reference, marginal, 1e-6);

  // derivatives w.r.t. the mean with more observations than latents
  constexpr stan::test::ad_tolerances tols{
      stan::test::ad_gradient_tols{1e-8, 1e-3}};
  auto f = [&](auto&& mean_arg) {
    return laplace_marginal_tol_poisson_log_lpmf(
        y, y_index, mean_arg, hessian_block_size,
        stan::math::test::squared_kernel_functor{},
        std::forward_as_tuple(x, alpha, rho),
        std::make_tuple(theta_0, tolerance, max_num_steps, solver,
                        max_steps_line_search, true),
        nullptr);
  };
  stan::test::expect_ad<true>(tols, f, mean);
}

}  // namespace
