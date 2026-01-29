#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <test/unit/math/laplace/aki_synth_data/x1.hpp>
#include <test/unit/math/laplace/motorcycle_gp/x_vec.hpp>
#include <test/unit/pretty_print_types.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <sstream>
#include <vector>

namespace {
struct poisson_log_likelihood2 {
  template <typename Theta>
  auto operator()(const Theta& theta, const std::vector<int>& delta_int,
                  std::ostream* pstream) const {
    return stan::math::poisson_log_lpmf(delta_int, theta);
  }
};

struct poisson_log_likelihood_tuple {
  template <typename Theta, typename Eta>
  auto operator()(const Theta& theta, const std::vector<int>& delta_int,
                  Eta&& eta, std::ostream* pstream) const {
    return stan::math::poisson_log_lpmf(delta_int, theta) + std::get<0>(eta)
           + std::get<1>(eta);
  }
};

struct poisson_log_likelihood_tuple_expanded {
  template <typename Theta, typename Eta, typename Eta1, typename Eta2>
  auto operator()(const Theta& theta, const std::vector<int>& delta_int,
                  Eta&& eta, Eta1&& eta1, Eta2&& eta2,
                  std::ostream* pstream) const {
    return stan::math::poisson_log_lpmf(delta_int, theta) + std::get<0>(eta)
           + std::get<1>(eta) + stan::math::sum(eta1) + stan::math::sum(eta2);
  }
};

struct poisson_re_log_ll {
  template <
      typename T0__, typename T2__,
      stan::require_all_t<
          stan::is_col_vector<T0__>, stan::is_vt_not_complex<T0__>,
          stan::is_col_vector<T2__>, stan::is_vt_not_complex<T2__>>* = nullptr>
  stan::return_type_t<stan::base_type_t<T0__>, stan::base_type_t<T2__>>
  operator()(const T0__& theta_arg__, const std::vector<int>& y,
             const T2__& mu_arg__, std::ostream* pstream__) const {
    const auto& theta = stan::math::to_ref(theta_arg__);
    const auto& mu = stan::math::to_ref(mu_arg__);
    return stan::math::poisson_log_lpmf<false>(y, stan::math::add(mu, theta));
  }
};

struct cov_fun {
  template <typename T0__,
            stan::require_all_t<stan::math::disjunction<
                stan::is_autodiff_scalar<T0__>,
                std::is_floating_point<std::decay_t<T0__>>>>* = nullptr>
  Eigen::Matrix<stan::return_type_t<T0__>, -1, -1> operator()(
      const T0__& sigma, const int& N, std::ostream* pstream__) const {
    return stan::math::diag_matrix(
        stan::math::rep_vector(stan::math::pow(sigma, 2), N));
  }
};

class laplace_types : public LaplaceAdTest {};

TEST(laplace, theta_0_as_expression_issue_3196) {
  // See https://github.com/stan-dev/math/issues/3196
  std::vector<int> y{1, 1, 1, 1, 1};
  Eigen::VectorXd beta{{1, 1, 1, 1, 1}};
  Eigen::VectorXd offset{{1, 1, 1, 1, 1}};
  double sigmaz = 1.0;
  double alpha = 1.0;
  int N = 5;
  Eigen::MatrixXd X{{{1, 1, 1, 1, 1},
                     {1, 1, 1, 1, 1},
                     {1, 1, 1, 1, 1},
                     {1, 1, 1, 1, 1},
                     {1, 1, 1, 1, 1}}};

  constexpr double tolerance = 1e-6;
  constexpr int max_num_steps = 1000;
  constexpr int hessian_block_size = 1;
  constexpr int solver_num = 1;
  constexpr int max_steps_line_search = 10;

  EXPECT_NO_THROW(stan::math::laplace_marginal_tol<false>(
      poisson_re_log_ll(),
      std::tuple<const std::vector<int>&, Eigen::Matrix<double, -1, 1>>(
          y, stan::math::add(stan::math::add(offset, alpha),
                             stan::math::multiply(X, beta))),
      cov_fun(), std::tuple<double, int>(sigmaz, N),
      std::tuple{stan::math::rep_vector(0.0, N), tolerance, max_num_steps,
                 hessian_block_size, solver_num, max_steps_line_search, 1},
      nullptr));
  auto arena_init = stan::math::to_arena(stan::math::rep_vector(0.0, N));
  EXPECT_NO_THROW(stan::math::laplace_marginal_tol<false>(
      poisson_re_log_ll(),
      std::tuple<const std::vector<int>&, Eigen::Matrix<double, -1, 1>>(
          y, stan::math::add(stan::math::add(offset, alpha),
                             stan::math::multiply(X, beta))),
      cov_fun(), std::tuple(sigmaz, N),
      std::make_tuple(arena_init, tolerance, max_num_steps, hessian_block_size,
                      solver_num, max_steps_line_search, true),
      nullptr));
}

TEST_P(laplace_types, poisson_log_phi_dim_2_tuple_extended) {
  using stan::math::laplace_marginal;
  using stan::math::laplace_marginal_tol;
  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;
  // logger->current_test_name_ = "poisson_log_phi_dim_2";
  constexpr int dim_phi = 2;
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_dbl(dim_phi);
  phi_dbl << 1.6, 0.45;

  int dim_theta = 2;
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
  const auto [solver_num, hessian_block_size, max_steps_line_search]
      = GetParam();
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);
  LAPLACE_SKIP_ZERO_STEPS(max_steps_line_search);

  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 1000;
  using stan::is_var_v;
  using stan::scalar_type_t;
  using stan::math::test::laplace_issue;
  constexpr stan::test::ad_tolerances tols{
      stan::test::ad_gradient_tols{1e-8, 1e-2}};
  //  stan::test::ad_tolerances tols;
  //  tols.gradient_grad_ = 1e-1;
  auto f_ll = [&](auto&& eta1, auto&& eta2, auto&& eta3) {
    try {
      auto eta1_tuple = std::make_tuple(eta1(0), eta1(1));
      return laplace_marginal_tol<false>(
          poisson_log_likelihood_tuple_expanded{},
          std::forward_as_tuple(sums, eta1_tuple, eta2, eta3),
          stan::math::test::squared_kernel_functor{},
          std::forward_as_tuple(x, std::make_tuple(phi_dbl(0), phi_dbl(1))),
          std::make_tuple(theta_0, tolerance, max_num_steps, hessian_block_size,
                          solver_num, max_steps_line_search, true),
          &output_stream);
    } catch (const std::exception& e) {
      std::stringstream fail_msg;
      using stan::math::test::test_type_name;
      fail_msg << "Exception thrown with eta1("
               << test_type_name<decltype(eta1)>() << ")=" << eta1[0]
               << ", eta2(" << test_type_name<decltype(eta2)>()
               << ")=" << eta2[0] << ", eta3("
               << test_type_name<decltype(eta3)>() << ")=" << eta3[0] << ". ";
      ADD_FAILURE() << fail_msg.str() << "\n Error message: " << e.what();
      throw;
    }
  };
  Eigen::VectorXd test1(phi_dbl);
  std::vector<double> test2 = {1.0, 1.0};
  stan::test::expect_ad<true>(tols, f_ll, phi_dbl, test1, test2);
}

TEST_P(laplace_types, poisson_log_phi_dim_2_tuple) {
  using stan::math::laplace_marginal;
  using stan::math::laplace_marginal_tol;
  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;
  // logger->current_test_name_ = "poisson_log_phi_dim_2";
  constexpr int dim_phi = 2;
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_dbl(dim_phi);
  phi_dbl << 1.6, 0.45;

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
  const auto [solver_num, hessian_block_size, max_steps_line_search]
      = GetParam();
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);
  LAPLACE_SKIP_ZERO_STEPS(max_steps_line_search);

  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 1000;
  using stan::is_var_v;
  using stan::scalar_type_t;
  using stan::math::test::laplace_issue;
  constexpr stan::test::ad_tolerances tols{
      stan::test::ad_gradient_tols{1e-8, 1e-3}};
  //  stan::test::ad_tolerances tols;
  //  tols.gradient_grad_ = 1e-1;
  auto f_covar = [&](auto&& x_v, auto&& alpha, auto&& rho) {
    try {
      return laplace_marginal_tol<false>(
          poisson_log_likelihood2{}, std::forward_as_tuple(sums),
          stan::math::test::squared_kernel_functor{},
          std::forward_as_tuple(x_v, std::make_tuple(alpha, rho)),
          std::make_tuple(theta_0, tolerance, max_num_steps, hessian_block_size,
                          solver_num, max_steps_line_search, true),
          &output_stream);
    } catch (const std::exception& e) {
      std::stringstream fail_msg;
      using stan::math::test::test_type_name;
      fail_msg << "Exception thrown with x(" << test_type_name<decltype(x_v)>()
               << "), alpha(" << test_type_name<decltype(alpha)>()
               << ")=" << alpha << ", rho(" << test_type_name<decltype(rho)>()
               << ")=" << rho << ". ";
      ADD_FAILURE() << fail_msg.str() << "\n Error message: " << e.what();
      throw;
    }
  };
  stan::test::expect_ad<true>(tols, f_covar, x, phi_dbl[0], phi_dbl[1]);
  auto f_ll = [&](auto&& alpha_rho, auto&& eta1, auto&& eta2) {
    try {
      return laplace_marginal_tol<false>(
          poisson_log_likelihood_tuple{},
          std::forward_as_tuple(sums, std::make_tuple(eta1, eta2)),
          stan::math::test::squared_kernel_functor{},
          std::forward_as_tuple(x, std::make_tuple(alpha_rho(0), alpha_rho(1))),
          std::make_tuple(theta_0, tolerance, max_num_steps, hessian_block_size,
                          solver_num, max_steps_line_search, true),
          &output_stream);
    } catch (const std::exception& e) {
      std::stringstream fail_msg;
      using stan::math::test::test_type_name;
      fail_msg << "Exception thrown with alpha_rho("
               << test_type_name<decltype(alpha_rho)>() << ")=" << alpha_rho
               << ", eta1(" << test_type_name<decltype(eta1)>() << ")=" << eta1
               << ", eta2(" << test_type_name<decltype(eta2)>() << ")=" << eta2
               << ". ";
      ADD_FAILURE() << fail_msg.str() << "\n Error message: " << e.what();
      throw;
    }
  };
  constexpr auto test1 = 1.0;
  constexpr auto test2 = 1.0;
  stan::test::expect_ad<true>(tols, f_ll, phi_dbl, test1, test2);
}

struct poisson_log_likelihood_array_tuple {
  template <typename Theta, typename Eta>
  auto operator()(const Theta& theta, const std::vector<int>& delta_int,
                  Eta&& eta, std::ostream* pstream) const {
    return stan::math::poisson_log_lpmf(delta_int, theta) + std::get<0>(eta[0])
           + std::get<1>(eta[0]);
  }
};

TEST_P(laplace_types, poisson_log_phi_dim_2_array_tuple) {
  using stan::math::laplace_marginal;
  using stan::math::laplace_marginal_tol;
  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;
  // logger->current_test_name_ = "poisson_log_phi_dim_2";
  constexpr int dim_phi = 2;
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_dbl(dim_phi);
  phi_dbl << 1.6, 0.45;

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
  const auto [solver_num, hessian_block_size, max_steps_line_search]
      = GetParam();
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);
  LAPLACE_SKIP_ZERO_STEPS(max_steps_line_search);

  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 1000;
  using stan::is_var_v;
  using stan::scalar_type_t;
  using stan::math::test::laplace_issue;
  constexpr stan::test::ad_tolerances tols{
      stan::test::ad_gradient_tols{1e-8, 1e-3}};
  //  stan::test::ad_tolerances tols;
  // tols.gradient_grad_ = 1e-1;
  auto f_ll = [&](auto&& alpha_rho, auto&& eta1, auto&& eta2) {
    try {
      std::vector<std::tuple<std::decay_t<decltype(eta1)>,
                             std::decay_t<decltype(eta2)>>>
          eta_tuple;
      eta_tuple.push_back(std::make_tuple(eta1, eta2));
      using alpha_scalar = stan::scalar_type_t<decltype(alpha_rho)>;
      std::vector<std::tuple<alpha_scalar, alpha_scalar>> alpha_tuple;
      alpha_tuple.push_back(std::make_tuple(alpha_rho(0), alpha_rho(1)));
      return laplace_marginal_tol<false>(
          poisson_log_likelihood_array_tuple{},
          std::forward_as_tuple(sums, eta_tuple),
          stan::math::test::squared_kernel_functor{},
          std::forward_as_tuple(x, alpha_tuple),
          std::make_tuple(theta_0, tolerance, max_num_steps, hessian_block_size,
                          solver_num, max_steps_line_search, true),
          &output_stream);
    } catch (const std::exception& e) {
      std::stringstream fail_msg;
      using stan::math::test::test_type_name;
      fail_msg << "Exception thrown with alpha_rho("
               << test_type_name<decltype(alpha_rho)>() << ")=" << alpha_rho
               << ", eta1(" << test_type_name<decltype(eta1)>() << ")=" << eta1
               << ", eta2(" << test_type_name<decltype(eta2)>() << ")=" << eta2
               << ". ";
      ADD_FAILURE() << fail_msg.str() << "\n Error message: " << e.what();
      throw;
    }
  };
  constexpr auto test1 = 1.0;
  constexpr auto test2 = 1.0;
  stan::test::expect_ad<true>(tols, f_ll, phi_dbl, test1, test2);
}

LAPLACE_INSTANTIATE_TEST_SUITE_P(laplace_types);

}  // namespace
