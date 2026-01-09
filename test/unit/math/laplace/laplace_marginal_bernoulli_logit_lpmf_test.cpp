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
  const auto [solver_num, hessian_block_size, max_steps_line_search]
      = GetParam();
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);
  LAPLACE_SKIP_ZERO_STEPS(max_steps_line_search);

  auto x1 = stan::test::laplace::x1;
  auto x2 = stan::test::laplace::x2;
  auto y = stan::test::laplace::y;
  std::vector<Eigen::VectorXd> x(dim_theta);
  for (int i = 0; i < dim_theta; i++) {
    x[i] = Eigen::VectorXd{{x1[i], x2[i]}};
  }
  std::vector<int> n_samples = stan::math::rep_array(1, dim_theta);
  Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);
  Eigen::VectorXd mean = Eigen::VectorXd::Zero(dim_theta);
  std::vector<double> delta;
  std::vector<int> delta_int;
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_dbl{{1.6, 1}};
  using stan::math::test::sqr_exp_kernel_functor;
  double target = laplace_marginal_bernoulli_logit_lpmf(
      y, n_samples, 0, sqr_exp_kernel_functor{},
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
          y, n_samples, mean, sqr_exp_kernel_functor{},
          std::forward_as_tuple(x, alpha, rho),
          std::make_tuple(theta_0, tolerance, max_num_steps, hessian_block_size,
                          solver_num, max_steps_line_search, true),
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

}  // namespace
