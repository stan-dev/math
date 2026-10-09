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
  const auto test_params = GetParam();
  const auto solver_num = std::get<0>(test_params);
  const auto hessian_block_size = std::get<1>(test_params);
  const auto max_steps_line_search = std::get<2>(test_params);
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);

  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 1000;
  constexpr stan::test::ad_tolerances tols{
      stan::test::ad_gradient_tols{1e-8, 1e-2}};
  auto f = [&](auto&& alpha, auto&& rho, auto&& eta) {
    try {
      return laplace_marginal_tol_neg_binomial_2_log_lpmf(
          y, y_index, eta, 0, hessian_block_size,
          stan::math::test::squared_kernel_functor{},
          std::forward_as_tuple(x, alpha, rho),
          std::make_tuple(theta_0, tolerance, max_num_steps, solver_num,
                          max_steps_line_search, true),
          &output_stream);
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
  const auto test_params = GetParam();
  const auto solver_num = std::get<0>(test_params);
  const auto hessian_block_size = std::get<1>(test_params);
  const auto max_steps_line_search = std::get<2>(test_params);
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);
  constexpr double eta = 1;

  // TODO(charlesm93): get benchmark from GPStuff or another software.
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 1000;
  auto smoke = [&](auto&& alpha, auto&& rho, auto&& eta_arg) {
    return laplace_marginal_tol_neg_binomial_2_log_lpmf(
        y, y_index, eta_arg, mean, hessian_block_size,
        stan::math::test::sqr_exp_kernel_functor{},
        std::forward_as_tuple(x, alpha, rho),
        std::make_tuple(theta_0, tolerance, max_num_steps, solver_num,
                        max_steps_line_search, true),
        &output_stream);
  };
  smoke(phi_dbl[0], phi_dbl[1], eta);
  auto f = [&](auto&& alpha, auto&& rho, auto&& eta_arg) {
    try {
      return laplace_marginal_tol_neg_binomial_2_log_lpmf(
          y, y_index, eta_arg, mean, hessian_block_size,
          stan::math::test::sqr_exp_kernel_functor{},
          std::forward_as_tuple(x, alpha, rho),
          std::make_tuple(theta_0, tolerance, max_num_steps, solver_num,
                          max_steps_line_search, true),
          &output_stream);
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

// Reference likelihood computed per observation with the prim lpmf,
// indexing the dispersion by group when it is a vector.
struct neg_binomial_2_log_obs_likelihood {
  template <typename Theta, typename Eta, typename Mean>
  auto operator()(const Theta& theta, const Eta& eta, const std::vector<int>& y,
                  const std::vector<int>& y_index, const Mean& mean,
                  std::ostream* /*pstream*/) const {
    Eigen::Matrix<stan::return_type_t<Theta, Mean>, Eigen::Dynamic, 1>
        theta_obs(y.size());
    Eigen::Matrix<stan::scalar_type_t<Eta>, Eigen::Dynamic, 1> eta_obs(
        y.size());
    for (size_t i = 0; i < y.size(); ++i) {
      theta_obs(i) = theta(y_index[i] - 1) + mean(y_index[i] - 1);
      if constexpr (stan::is_stan_scalar<Eta>::value) {
        eta_obs(i) = eta;
      } else {
        eta_obs(i) = eta(y_index[i] - 1);
      }
    }
    return stan::math::neg_binomial_2_log_lpmf(y, theta_obs, eta_obs);
  }
};

class laplace_marginal_neg_binomial_log_lpmf_grouped : public ::testing::Test {
 protected:
  const std::vector<int> y{1, 0, 5, 2, 7, 0};
  const std::vector<int> y_index{1, 2, 3, 1, 2, 1};
  const Eigen::VectorXd mean{{0.3, -0.2, 0.1}};
  const std::vector<Eigen::VectorXd> x{
      Eigen::VectorXd{{0.05100797, 0.16086164}},
      Eigen::VectorXd{{-0.59823393, 0.98701425}},
      Eigen::VectorXd{{0.31296868, -0.68926772}}};
  const Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(3);
  static constexpr double alpha = 1.6;
  static constexpr double rho = 0.45;
  static constexpr double tolerance = 1e-12;
  static constexpr int max_num_steps = 1000;
  static constexpr int hessian_block_size = 1;
  static constexpr int solver = 1;
  static constexpr int max_steps_line_search = 0;

  template <typename Eta>
  auto marginal(const Eta& eta) {
    return stan::math::laplace_marginal_tol_neg_binomial_2_log_lpmf(
        y, y_index, eta, mean, hessian_block_size,
        stan::math::test::squared_kernel_functor{},
        std::forward_as_tuple(x, alpha, rho),
        std::make_tuple(theta_0, tolerance, max_num_steps, solver,
                        max_steps_line_search, true),
        nullptr);
  }
  template <typename Eta>
  auto reference(const Eta& eta) {
    return stan::math::laplace_marginal_tol<false>(
        neg_binomial_2_log_obs_likelihood{},
        std::forward_as_tuple(eta, y, y_index, mean), hessian_block_size,
        stan::math::test::squared_kernel_functor{},
        std::forward_as_tuple(x, alpha, rho),
        std::make_tuple(theta_0, tolerance, max_num_steps, solver,
                        max_steps_line_search, true),
        nullptr);
  }
};

// Scalar dispersion with more observations than latent variables.
TEST_F(laplace_marginal_neg_binomial_log_lpmf_grouped, scalar_eta) {
  constexpr double eta = 1.5;
  EXPECT_NEAR(reference(eta), marginal(eta), 1e-6);
}

// A vector dispersion (one entry per group) must work and match the
// per-observation reference. This is what stanc generates for the Stan
// signature, which requires `vector eta`.
TEST_F(laplace_marginal_neg_binomial_log_lpmf_grouped, vector_eta) {
  const Eigen::VectorXd eta{{1.5, 2.5, 0.75}};
  EXPECT_NEAR(reference(eta), marginal(eta), 1e-6);

  // a constant vector dispersion agrees with the scalar version
  const Eigen::VectorXd eta_rep = Eigen::VectorXd::Constant(3, 1.5);
  EXPECT_NEAR(marginal(1.5), marginal(eta_rep), 1e-8);

  // dispersion vector must have one entry per latent variable
  const Eigen::VectorXd eta_bad = Eigen::VectorXd::Constant(2, 1.5);
  EXPECT_THROW(marginal(eta_bad), std::invalid_argument);
}

// eta and mean must remain differentiable (autodiff arguments).
TEST_F(laplace_marginal_neg_binomial_log_lpmf_grouped, vector_eta_ad) {
  const Eigen::VectorXd eta{{1.5, 2.5, 0.75}};
  constexpr stan::test::ad_tolerances tols{
      stan::test::ad_gradient_tols{1e-8, 1e-3}};
  auto f = [&](auto&& eta_arg, auto&& mean_arg) {
    return stan::math::laplace_marginal_tol_neg_binomial_2_log_lpmf(
        y, y_index, eta_arg, mean_arg, hessian_block_size,
        stan::math::test::squared_kernel_functor{},
        std::forward_as_tuple(x, alpha, rho),
        std::make_tuple(theta_0, tolerance, max_num_steps, solver,
                        max_steps_line_search, true),
        nullptr);
  };
  stan::test::expect_ad<true>(tols, f, eta, mean);
}

}  // namespace
