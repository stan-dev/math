#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/laplace/aki_synth_data/x1.hpp>
#include <test/unit/pretty_print_types.hpp>

#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <sstream>
#include <vector>

namespace {

class laplace_marginal_bernoulli_logit_lpmf : public LaplaceAdTest {};

TEST_P(laplace_marginal_bernoulli_logit_lpmf, phi_dim500) {
  using stan::math::laplace_marginal_bernoulli_logit_lpmf;
  using stan::math::laplace_marginal_tol_bernoulli_logit_lpmf;
  using stan::math::to_vector;
  using stan::math::var;
  using stan::math::test::flag_test;
  constexpr int dim_theta = 500;
  const auto test_params = GetParam();
  const auto solver_num = std::get<0>(test_params);
  const auto hessian_block_size = std::get<1>(test_params);
  const auto max_steps_line_search = std::get<2>(test_params);
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);
  // LAPLACE_SKIP_ZERO_STEPS(max_steps_line_search);

  auto x1 = stan::test::laplace::x1;
  auto x2 = stan::test::laplace::x2;
  auto y = stan::test::laplace::y;
  std::vector<Eigen::VectorXd> x(dim_theta);
  for (int i = 0; i < dim_theta; i++) {
    x[i] = Eigen::VectorXd{{x1[i], x2[i]}};
  }
  std::vector<int> y_index;
  y_index.reserve(dim_theta);
  for (int i = 1; i <= dim_theta; i++) {
    y_index.push_back(i);
  }
  Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);
  Eigen::VectorXd mean = Eigen::VectorXd::Zero(dim_theta);
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_dbl{{1.6, 1}};
  using stan::math::test::sqr_exp_kernel_functor;
  double target = laplace_marginal_bernoulli_logit_lpmf(
      y, y_index, 0, hessian_block_size, sqr_exp_kernel_functor{},
      std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1)), nullptr);
  // Benchmark against gpstuff.
  constexpr double tol = 8e-4;
  EXPECT_NEAR(-195.368, target, tol);
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 1000;
  constexpr stan::test::ad_tolerances tols{
      stan::test::ad_gradient_tols{1e-8, 1e-3}};
  auto f = [&](auto&& alpha, auto&& rho) {
    try {
      return laplace_marginal_tol_bernoulli_logit_lpmf(
          y, y_index, mean, hessian_block_size, sqr_exp_kernel_functor{},
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
  stan::test::expect_ad<true>(tols, f, phi_dbl[0], phi_dbl[1]);
}

LAPLACE_INSTANTIATE_TEST_SUITE_P(laplace_marginal_bernoulli_logit_lpmf);

// Reference likelihood computed per observation, without the grouped
// sufficient-statistics shortcut used by bernoulli_logit_likelihood.
struct bernoulli_logit_obs_likelihood {
  template <typename Theta, typename Mean>
  auto operator()(const Theta& theta, const std::vector<int>& y,
                  const std::vector<int>& y_index, const Mean& mean,
                  std::ostream* /*pstream*/) const {
    Eigen::Matrix<stan::return_type_t<Theta, Mean>, Eigen::Dynamic, 1>
        theta_obs(y.size());
    for (size_t i = 0; i < y.size(); ++i) {
      theta_obs(i) = theta(y_index[i] - 1) + mean(y_index[i] - 1);
    }
    return stan::math::bernoulli_logit_lpmf(y, theta_obs);
  }
};

// More observations than latent variables: groups with several observations
// must aggregate all of them, not just the first theta.size() entries.
TEST(laplace_marginal_bernoulli_logit_lpmf_grouped, multiple_obs_per_group) {
  using stan::math::laplace_marginal_tol;
  using stan::math::laplace_marginal_tol_bernoulli_logit_lpmf;

  constexpr int dim_theta = 3;
  const std::vector<int> y{1, 0, 1, 1, 1, 0};
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

  const double marginal = laplace_marginal_tol_bernoulli_logit_lpmf(
      y, y_index, mean, hessian_block_size,
      stan::math::test::squared_kernel_functor{},
      std::forward_as_tuple(x, alpha, rho),
      std::make_tuple(theta_0, tolerance, max_num_steps, solver,
                      max_steps_line_search, true),
      nullptr);
  const double reference = laplace_marginal_tol<false>(
      bernoulli_logit_obs_likelihood{}, std::forward_as_tuple(y, y_index, mean),
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
    return laplace_marginal_tol_bernoulli_logit_lpmf(
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
