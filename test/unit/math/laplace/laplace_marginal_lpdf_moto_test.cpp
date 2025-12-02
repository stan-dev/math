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

struct normal_likelihood {
  template <typename Theta, typename YVec>
  auto operator()(const Theta& theta, const YVec& y, const int delta_int,
                  std::ostream* pstream) const {
    int n_obs = delta_int;
    Eigen::Matrix<stan::return_type_t<Theta>, -1, 1> mu(n_obs);
    Eigen::Matrix<stan::return_type_t<Theta>, -1, 1> sigma(n_obs);
    stan::return_type_t<Theta> lp = 0;
    for (Eigen::Index i = 0; i < n_obs; i++) {
      mu(i) = theta(2 * i);
      // TODO(Charles): Theta can be a large negative value so sigma can be 0
      sigma(i) = stan::math::lb_constrain<true>(
          stan::math::multiply(0.5, theta(2 * i + 1)), 1e-14, lp);
    }
    try {
      return stan::math::normal_lpdf(y, mu, sigma) + lp;
    } catch (const std::domain_error& e) {
      std::cout << "Error in normal_lpdf: " << e.what() << std::endl;
      std::cout << "theta: \n" << theta.transpose() << std::endl;
      std::cout << "y: \n" << y.transpose() << std::endl;
      std::cout << "mu: \n" << mu.transpose() << std::endl;
      std::cout << "sigma: \n" << sigma.transpose() << std::endl;
      return stan::math::normal_lpdf(y, mu, sigma) + lp;
    }
  }
};

struct covariance_motorcycle_functor {
  template <typename TX, typename LengthF, typename LengthG, typename SigmaF,
            typename SigmaG>
  auto operator()(const TX& x, const LengthF& length_scale_f,
                  const LengthG& length_scale_g, const SigmaF& sigma_f,
                  const SigmaG& sigma_g, const int n_obs,
                  std::ostream* msgs = nullptr) const {
    using Eigen::Matrix;
    using stan::math::gp_exp_quad_cov;
    using scalar_t = stan::return_type_t<LengthF, LengthG, SigmaF, SigmaG>;

    constexpr double jitter = 1e-8;
    Matrix<scalar_t, -1, -1> kernel_f
        = gp_exp_quad_cov(x, sigma_f, length_scale_f);
    Matrix<scalar_t, -1, -1> kernel_g
        = gp_exp_quad_cov(x, sigma_g, length_scale_g);

    Matrix<scalar_t, -1, -1> kernel_all
        = Eigen::MatrixXd::Zero(2 * n_obs, 2 * n_obs);
    for (Eigen::Index i = 0; i < n_obs; i++) {
      for (Eigen::Index j = 0; j <= i; j++) {
        kernel_all(2 * i, 2 * j) = kernel_f(i, j);
        kernel_all(2 * i + 1, 2 * j + 1) = kernel_g(i, j);
        if (i != j) {
          kernel_all(2 * j, 2 * i) = kernel_all(2 * i, 2 * j);
          kernel_all(2 * j + 1, 2 * i + 1) = kernel_all(2 * i + 1, 2 * j + 1);
        }
      }
    }
    for (Eigen::Index i = 0; i < 2 * n_obs; i++) {
      kernel_all(i, i) += jitter;
    }
    return kernel_all;
  }
};

class laplace_motorcyle_gp_test : public LaplaceAdTest {
 protected:
  void SetUp() override {
    using stan::math::gp_exp_quad_cov;
    using stan::math::value_of;
    Eigen::MatrixXd K_plus_I
        = gp_exp_quad_cov(x, value_of(sigma_f), value_of(length_scale_f))
          + Eigen::MatrixXd::Identity(n_obs, n_obs);
    Eigen::VectorXd mu_hat = K_plus_I.colPivHouseholderQr().solve(y);
    // Remark: finds optimal point with or without informed initial guess.
    // Better θ0: μ at GP posterior mean; g at a stable constant σ from
    // residuals
    for (int i = 0; i < n_obs; ++i) {
      theta0(2 * i) = mu_hat(i);
    }
    // After computing mu_hat as you already do
    Eigen::VectorXd r = (y - mu_hat).cwiseAbs();

    // Optional tiny smoothing to avoid zeros / spikes (window radius = 2)
    Eigen::VectorXd r_smooth = r;
    for (int i = 0; i < n_obs; ++i) {
      double acc = 0.0;
      int cnt = 0;
      for (int j = std::max(0, i - 2); j <= std::min(n_obs - 1, i + 2); ++j) {
        acc += r(j);
        ++cnt;
      }
      r_smooth(i) = acc / cnt;
    }

    // Baseline scale for clamping
    double s0 = std::max(1e-3, std::sqrt(r.array().square().mean()));

    // Choose sigma_i0 close to |residual|, and clamp to a sane band
    for (int i = 0; i < n_obs; ++i) {
      double si = std::min(2.0 * s0, std::max(0.5 * s0, r_smooth(i)));
      // sigma = lb + exp(0.5 * theta)  =>  theta = 2 * log(sigma - lb)
      theta0(2 * i + 1) = 2.0 * std::log(std::max(si - 1e-14, 1e-12));
    }
  }

  static constexpr int n_obs{133};
  static constexpr int dim_phi{4};
  std::vector<double> x{stan::test::laplace::moto::x};
  Eigen::VectorXd y{stan::test::laplace::moto::y};

  static constexpr double length_scale_f{0.3};
  static constexpr double length_scale_g{0.5};
  static constexpr double sigma_f{0.25};
  static constexpr double sigma_g{0.25};
  Eigen::VectorXd theta0{Eigen::VectorXd::Zero(2 * n_obs)};
  Eigen::Matrix<double, -1, 1> eta{{1.0}};
  static constexpr double sigma_global{1.0};
  Eigen::VectorXd eta_dbl{{sigma_global}};
  static constexpr int solver{2};
  static constexpr double eps{1e-7};
  Eigen::VectorXd phi_dbl{{length_scale_f, length_scale_g, sigma_f, sigma_g}};
};

TEST_P(laplace_motorcyle_gp_test, gp_motorcycle_val) {
  // logger->current_test_name_ = "gp_motorcycle";
  using stan::math::laplace_marginal_tol;
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 1000;
  const auto [solver_num, hessian_block_size, max_steps_line_search]
      = GetParam();
  constexpr int dim_theta = 2 * n_obs;
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);
  LAPLACE_SKIP_ZERO_STEPS(max_steps_line_search);

  double target = laplace_marginal_tol<false>(
      normal_likelihood{}, std::forward_as_tuple(y, n_obs),
      covariance_motorcycle_functor{},
      std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1), phi_dbl(2), phi_dbl(3),
                            n_obs),
      theta0, tolerance, max_num_steps, hessian_block_size, solver_num,
      max_steps_line_search, nullptr);
}

TEST_P(laplace_motorcyle_gp_test, gp_motorcycle_ad) {
  using stan::math::gp_exp_quad_cov;
  using stan::math::value_of;
  // logger->current_test_name_ = "gp_motorcycle";
  using stan::math::laplace_marginal_tol;

  // TODO(Steve): benchmark this result against GPStuff.
  constexpr double tolerance = 1e-8;
  constexpr int max_num_steps = 1000;
  auto phi_0 = phi_dbl(0);
  auto phi_1 = phi_dbl(1);
  Eigen::VectorXd phi_rest = phi_dbl.tail(2);
  Eigen::VectorXd phi_01{{phi_0, phi_1}};
  const auto [solver_num, hessian_block_size, max_steps_line_search]
      = GetParam();
  constexpr int dim_theta = 2 * n_obs;
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);
  LAPLACE_SKIP_ZERO_STEPS(max_steps_line_search);
  constexpr stan::test::ad_tolerances tols{
      stan::test::ad_gradient_tols{1e-8, 1e-1}};
  auto f = [&](auto&& phi_01_v, auto&& phi_rest_v) {
    try {
      return laplace_marginal_tol<false>(
          normal_likelihood{}, std::forward_as_tuple(y, n_obs),
          covariance_motorcycle_functor{},
          std::forward_as_tuple(x, phi_01_v(0), phi_01_v(1), phi_rest_v(0),
                                phi_rest_v(1), n_obs),
          theta0, tolerance, max_num_steps, hessian_block_size, solver_num,
          max_steps_line_search, nullptr);
    } catch (const std::exception& e) {
      std::stringstream fail_msg;
      using stan::math::test::test_type_name;
      fail_msg << "Exception thrown with phi_01_v("
               << test_type_name<decltype(phi_01_v)>() << ")=" << phi_01_v
               << ", phi_rest_v(" << test_type_name<decltype(phi_rest_v)>()
               << ")=" << phi_rest_v << ". ";
      ADD_FAILURE() << fail_msg.str() << "\n Error message: " << e.what();
      throw;
    }
  };
  try {
    stan::test::expect_ad<true>(tols, f, phi_01, phi_rest);
  } catch (const std::domain_error e) {
    ADD_FAILURE() << "Exception: " << e.what()
                  << "\n\tsolver_num: " << solver_num
                  << "\n\tmax_steps_line_search: " << max_steps_line_search
                  << "\n\thessian_block_size: " << hessian_block_size
                  << std::endl;
    stan::math::recover_memory();
  }
}

struct normal_likelihood2 {
  template <typename Theta, typename SigmaGlobal>
  auto operator()(const Theta& theta, const Eigen::VectorXd& y, const int n_obs,
                  const SigmaGlobal& sigma_global,
                  std::ostream* pstream) const {
    using stan::math::multiply;
    Eigen::Matrix<stan::return_type_t<Theta>, -1, 1> mu(n_obs);
    Eigen::Matrix<stan::return_type_t<Theta>, -1, 1> sigma(n_obs);
    stan::return_type_t<Theta> lp = 0;
    for (int i = 0; i < n_obs; i++) {
      mu(i) = theta(2 * i);
      sigma(i) = stan::math::lb_constrain<true>(multiply(0.5, theta(2 * i + 1)),
                                                0.0, lp);  // * sigma_global;
    }
    // return stan::math::normal_lpdf(y, mu, sigma);
    return stan::math::normal_lpdf(y, mu, stan::math::add(sigma_global, sigma))
           + lp;
  }
};

TEST_P(laplace_motorcyle_gp_test, gp_motorcycle2_val) {
  using stan::math::gp_exp_quad_cov;
  using stan::math::value_of;
  Eigen::MatrixXd K_plus_I
      = gp_exp_quad_cov(x, value_of(sigma_f), value_of(length_scale_f))
        + Eigen::MatrixXd::Identity(n_obs, n_obs);
  Eigen::VectorXd mu_hat = K_plus_I.colPivHouseholderQr().solve(y);
  // Remark: finds optimal point with or without informed initial guess.
  for (int i = 0; i < n_obs - 1; i++) {
    theta0(2 * i) = mu_hat(i);
    theta0(2 * i + 1) = -1.0;
  }
  using stan::math::laplace_marginal_tol;
  Eigen::VectorXd length_scale_vec = phi_dbl.head(2);
  Eigen::VectorXd sigma_vec = phi_dbl.tail(2);
  constexpr double tolerance = 1e-12;
  constexpr int max_num_steps = 300;
  const auto [solver_num, hessian_block_size, max_steps_line_search]
      = GetParam();
  constexpr int dim_theta = 2 * n_obs;
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);
  LAPLACE_SKIP_ZERO_STEPS(max_steps_line_search);
  double target = laplace_marginal_tol<false>(
      normal_likelihood2{}, std::forward_as_tuple(y, n_obs, sigma_global),
      covariance_motorcycle_functor{},
      std::forward_as_tuple(x, length_scale_f, length_scale_g, sigma_f, sigma_g,
                            n_obs),
      theta0, tolerance, max_num_steps, hessian_block_size, solver_num,
      max_steps_line_search, nullptr);
}

TEST_P(laplace_motorcyle_gp_test, gp_motorcycle2_ad) {
  using stan::math::gp_exp_quad_cov;
  using stan::math::laplace_marginal_tol;
  using stan::math::value_of;
  Eigen::MatrixXd K_plus_I
      = gp_exp_quad_cov(x, value_of(sigma_f), value_of(length_scale_f))
        + Eigen::MatrixXd::Identity(n_obs, n_obs);
  Eigen::VectorXd mu_hat = K_plus_I.colPivHouseholderQr().solve(y);
  // TODO(Charles): benchmark this result against GPStuff.
  constexpr double tolerance = 1e-8;
  constexpr int max_num_steps = 1000;
  Eigen::VectorXd length_scale_vec = phi_dbl.head(2);
  Eigen::VectorXd sigma_vec = phi_dbl.tail(2);
  const auto [solver_num, hessian_block_size, max_steps_line_search]
      = GetParam();
  constexpr int dim_theta = 2 * n_obs;
  LAPLACE_SKIP_IF_INVALID_TEST_COMBO(hessian_block_size, dim_theta);
  LAPLACE_SKIP_ZERO_STEPS(max_steps_line_search);
  constexpr stan::test::ad_tolerances tols{
      stan::test::ad_gradient_tols{1e-8, 1e-1}};
  auto f = [&](auto&& sigma_global_v, auto&& length_scale_v, auto&& sigma_v) {
    try {
      return laplace_marginal_tol<false>(
          normal_likelihood2{}, std::forward_as_tuple(y, n_obs, sigma_global_v),
          covariance_motorcycle_functor{},
          std::forward_as_tuple(x, length_scale_v(0), length_scale_v(1),
                                sigma_v(0), sigma_v(1), n_obs),
          theta0, tolerance, max_num_steps, hessian_block_size, solver_num,
          max_steps_line_search, nullptr);
    } catch (const std::exception& e) {
      std::stringstream fail_msg;
      using stan::math::test::test_type_name;
      fail_msg << "Exception thrown with sigma_global_v("
               << test_type_name<decltype(sigma_global_v)>()
               << ")=" << sigma_global_v << ", length_scale_v("
               << test_type_name<decltype(length_scale_v)>()
               << ")=" << length_scale_v << ", sigma_v("
               << test_type_name<decltype(sigma_v)>() << ")=" << sigma_v
               << ". ";
      ADD_FAILURE() << fail_msg.str() << "\n Error message: " << e.what();
      throw;
    }
  };
  stan::test::expect_ad<true>(tols, f, sigma_global, length_scale_vec,
                              sigma_vec);
}

LAPLACE_INSTANTIATE_TEST_SUITE_P(laplace_motorcyle_gp_test);

}  // namespace
