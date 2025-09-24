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
  int dim_phi = 2;
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

  double target = laplace_marginal<false>(
      poisson_log_likelihood2{}, std::forward_as_tuple(sums),
      stan::math::test::squared_kernel_functor{},
      std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1)), nullptr);

  // TODO(Charles): benchmark target against gpstuff.
  // Expected: -2.53056
  double tol = 1e-4;
  EXPECT_NEAR(-2.53056, value_of(target), tol);

  // Test with optional arguments
  {
    constexpr double tolerance = 1e-8;
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
  using stan::math::test::laplace_issue;
  constexpr std::array known_issues{laplace_issue{0, 0, 0}};
  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = 1e-1;
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

    double tol = 6e-4;
    // Benchmark from GPStuff.
    EXPECT_NEAR(-2866.88, value_of(marginal_density), tol);
  }
  constexpr double tolerance = 1e-8;
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
  int dim_theta = 500;
  int n_observations = 500;
  auto x1 = stan::test::laplace::x1;
  auto x2 = stan::test::laplace::x2;
  auto y = stan::test::laplace::y;

  int dim_x = 2;
  std::vector<Eigen::VectorXd> x(dim_theta);
  for (int i = 0; i < dim_theta; i++) {
    Eigen::VectorXd coordinate(dim_x);
    coordinate << x1[i], x2[i];
    x[i] = coordinate;
  }
  Eigen::VectorXd theta_0 = Eigen::VectorXd::Zero(dim_theta);
  Eigen::VectorXd delta_L;
  std::vector<double> delta;
  int dim_phi = 2;
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_dbl(dim_phi);
  phi_dbl << 1.6, 1;

  double target = laplace_marginal<false>(
      bernoulli_logit_likelihood{}, std::forward_as_tuple(y),
      stan::math::test::sqr_exp_kernel_functor{},
      std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1)), nullptr);

  double tol = 8e-5;
  // Benchmark against gpstuff.
  EXPECT_NEAR(-195.368, target, tol);
  // All fail for ad check with relative tolerance ~0.002
  constexpr double tolerance = 1e-8;
  constexpr int max_num_steps = 100;
  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = 1e-3;
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

    double jitter = 1e-12;
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

struct normal_likelihood {
  template <typename Theta, typename YVec>
  auto operator()(const Theta& theta, const YVec& y, const int delta_int,
                  std::ostream* pstream) const {
    int n_obs = delta_int;
    Eigen::Matrix<stan::return_type_t<Theta>, -1, 1> mu(n_obs);
    Eigen::Matrix<stan::return_type_t<Theta>, -1, 1> sigma(n_obs);
    for (Eigen::Index i = 0; i < n_obs; i++) {
      mu(i) = theta(2 * i);
      // TODO(Charles): Theta can be a large negative value so sigma can be 0
      sigma(i) = exp(0.5 * theta(2 * i + 1)) + 1e-12;
    }
    try {
      return stan::math::normal_lpdf(y, mu, sigma);
    } catch (const std::domain_error& e) {
      std::cout << "Error in normal_lpdf: " << e.what() << std::endl;
      std::cout << "theta: \n" << theta.transpose() << std::endl;
      std::cout << "y: \n" << y.transpose() << std::endl;
      std::cout << "mu: \n" << mu.transpose() << std::endl;
      std::cout << "sigma: \n" << sigma.transpose() << std::endl;
      return stan::math::normal_lpdf(y, mu, sigma);
    }
  }
};

class laplace_motorcyle_gp_test : public ::testing::Test {
 protected:
  void SetUp() override {
    using stan::math::gp_exp_quad_cov;
    using stan::math::value_of;
    Eigen::MatrixXd K_plus_I
        = gp_exp_quad_cov(x, value_of(sigma_f), value_of(length_scale_f))
          + Eigen::MatrixXd::Identity(n_obs, n_obs);
    Eigen::VectorXd mu_hat = K_plus_I.colPivHouseholderQr().solve(y);
    // Remark: finds optimal point with or without informed initial guess.
    for (int i = 0; i < n_obs; i++) {
      theta0(2 * i) = mu_hat(i);  // 0
      theta0(2 * i + 1) = 0.1;
    }
  }

  int n_obs{133};
  int dim_phi{4};
  std::vector<double> x{stan::test::laplace::moto::x};
  Eigen::VectorXd y{stan::test::laplace::moto::y};

  double length_scale_f{0.3};
  double length_scale_g{0.5};
  double sigma_f{0.25};
  double sigma_g{0.25};
  std::vector<int> delta_int{n_obs};
  Eigen::VectorXd theta0{Eigen::VectorXd::Zero(2 * n_obs)};
  Eigen::Matrix<double, -1, 1> eta{{1.0}};
  Eigen::VectorXd eta_dbl{{1.0}};
  int solver{2};
  double eps{1e-7};
  Eigen::VectorXd phi_dbl{{length_scale_f, length_scale_g, sigma_f, sigma_g}};
};

TEST_F(laplace_motorcyle_gp_test, gp_motorcycle) {
  // logger->current_test_name_ = "gp_motorcycle";
  using stan::math::laplace_marginal;
  using stan::math::laplace_marginal_tol;
  using stan::math::value_of;

  {
    constexpr double tolerance = 1e-08;
    constexpr int max_num_steps = 100;
    constexpr int hessian_block_size = 2;
    solver = 2;
    constexpr int do_line_search = 1;
    constexpr int max_steps_line_search = 10;

    double target = laplace_marginal_tol<false>(
        normal_likelihood{}, std::forward_as_tuple(y, delta_int[0]),
        covariance_motorcycle_functor{},
        std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1), phi_dbl(2), phi_dbl(3),
                              n_obs),
        theta0, tolerance, max_num_steps, hessian_block_size, solver,
        max_steps_line_search, nullptr);
  }

  // TODO(Steve): benchmark this result against GPStuff.
  constexpr double tolerance = 1e-6;
  constexpr int max_num_steps = 1000;
  auto phi_0 = phi_dbl(0);
  auto phi_1 = phi_dbl(1);
  Eigen::VectorXd phi_rest = phi_dbl.tail(2);
  Eigen::VectorXd phi_01{{phi_0, phi_1}};
  using stan::math::test::laplace_issue;
  using stan::math::test::LaplaceFailures;
  constexpr std::array known_issues{
      std::pair(laplace_issue{1, 0, 1}, LaplaceFailures::HessianFailure),
      std::pair(laplace_issue{1, 100, 1}, LaplaceFailures::HessianFailure),
      std::pair(laplace_issue{1, 200, 1}, LaplaceFailures::HessianFailure),
      std::pair(laplace_issue{1, 300, 1}, LaplaceFailures::HessianFailure),
      std::pair(laplace_issue{1, 400, 1}, LaplaceFailures::HessianFailure),
      std::pair(laplace_issue{1, 500, 1}, LaplaceFailures::HessianFailure),
      std::pair(laplace_issue{1, 0, 2}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 100, 2}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 200, 2}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 300, 2}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 400, 2}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 500, 2}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 0, 3}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 100, 3}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 200, 3}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 300, 3}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 400, 3}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 500, 3}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 0, 4}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 100, 4}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 200, 4}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 300, 4}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 400, 4}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{1, 500, 4}, LaplaceFailures::SqrtDNE),
      std::pair(laplace_issue{2, 0, 1}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{2, 100, 1}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{2, 200, 1}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{2, 300, 1}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{2, 400, 1}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{2, 500, 1}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{2, 0, 3}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{2, 100, 3}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{2, 200, 3}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{2, 300, 3}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{2, 400, 3}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{2, 500, 3}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{2, 0, 4}, LaplaceFailures::IterExceeded),
      std::pair(laplace_issue{2, 100, 4}, LaplaceFailures::IterExceeded),
      std::pair(laplace_issue{2, 200, 4}, LaplaceFailures::IterExceeded),
      std::pair(laplace_issue{2, 300, 4}, LaplaceFailures::IterExceeded),
      std::pair(laplace_issue{2, 400, 4}, LaplaceFailures::IterExceeded),
      std::pair(laplace_issue{2, 500, 4}, LaplaceFailures::IterExceeded),
      std::pair(laplace_issue{3, 0, 1}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{3, 100, 1}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{3, 200, 1}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{3, 300, 1}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{3, 400, 1}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{3, 500, 1}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{3, 0, 3}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{3, 100, 3}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{3, 200, 3}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{3, 300, 3}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{3, 400, 3}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{3, 500, 3}, LaplaceFailures::NaNTheta),
      std::pair(laplace_issue{3, 0, 4}, LaplaceFailures::IterExceeded),
      std::pair(laplace_issue{3, 100, 4}, LaplaceFailures::IterExceeded),
      std::pair(laplace_issue{3, 200, 4}, LaplaceFailures::IterExceeded),
      std::pair(laplace_issue{3, 300, 4}, LaplaceFailures::IterExceeded),
      std::pair(laplace_issue{3, 400, 4}, LaplaceFailures::IterExceeded),
      std::pair(laplace_issue{3, 500, 4}, LaplaceFailures::IterExceeded)};

  /**
   * Note: This test is designed to check the error behavior
   *  of the laplace_marginal_tol function. We do not force
   *  a function to fail because some of these errors can be machine
   *  specific. So for cases we know there can be a test failure for a
   *  machine we call the function in a try block. if it *does* fail,
   *  we expect it to be the associated error found in the known_issues array.
   *  If we have not seen this parameter combination fail before, we run the
   *  standard AD testing procedure.
   */
  for (int solver_num = 1; solver_num < 4; solver_num++) {
    for (int max_steps_line_search = 0; max_steps_line_search <= 20;
         max_steps_line_search += 10) {
      for (int hessian_block_size = 1; hessian_block_size < 3;
           hessian_block_size++) {
        // logger->update_laplace_info(solver_num, hessian_block_size,
        // max_steps_line_search);
        if (theta0.size() % hessian_block_size != 0) {
          std::cerr << "[          ] [ INFO ]"
                    << " Skipping test for hessian of size " << theta0.size()
                    << " with hessian block size of " << hessian_block_size
                    << std::endl;
          continue;
        }
        auto f = [&](auto&& y_v, auto&& phi_01_v, auto&& phi_rest_v) {
          return laplace_marginal_tol<false>(
              normal_likelihood{}, std::forward_as_tuple(y_v, delta_int[0]),
              covariance_motorcycle_functor{},
              std::forward_as_tuple(x, phi_01_v(0), phi_01_v(0), phi_rest_v(0),
                                    phi_rest_v(1), n_obs),
              theta0, tolerance, max_num_steps, hessian_block_size, solver_num,
              max_steps_line_search, nullptr);
        };
        stan::test::ad_tolerances tols;
        tols.gradient_grad_ = 1e-1;
        using stan::math::test::flag_test;
        auto flag_val = flag_test(known_issues, solver_num,
                                  max_steps_line_search, hessian_block_size);
        if (flag_val != LaplaceFailures::None) {
          try {
            auto ret = f(y, phi_01, phi_rest);
          } catch (const std::domain_error& e) {
            using stan::math::test::err_to_laplace_failure;
            LaplaceFailures err_val = err_to_laplace_failure(e);
            EXPECT_EQ(err_val, flag_val)
                << "Error: " << e.what()
                << "\n\terr_val: " << to_string(err_val)
                << "\n\tflag_val: " << to_string(flag_val)
                << "\n\tsolver_num: " << solver_num
                << "\n\tmax_steps_line_search: " << max_steps_line_search
                << "\n\thessian_block_size: " << hessian_block_size;
          }
          stan::math::recover_memory();
        } else {
          try {
            stan::test::expect_ad<true>(tols, f, y, phi_01, phi_rest);
          } catch (const std::domain_error e) {
            ADD_FAILURE() << "Exception: " << e.what()
                          << "\n\tsolver_num: " << solver_num
                          << "\n\tmax_steps_line_search: "
                          << max_steps_line_search
                          << "\n\thessian_block_size: " << hessian_block_size
                          << std::endl;
            stan::math::recover_memory();
          }
        }
      }
    }
  }
}

struct normal_likelihood2 {
  template <typename Theta, typename Eta>
  auto operator()(const Theta& theta, const Eigen::VectorXd& y,
                  const std::vector<int>& delta_int, const Eta& eta,
                  std::ostream* pstream) const {
    using stan::math::multiply;
    int n_obs = delta_int[0];
    Eigen::Matrix<stan::return_type_t<Theta>, -1, 1> mu(n_obs);
    Eigen::Matrix<stan::return_type_t<Theta>, -1, 1> sigma(n_obs);
    auto sigma_global = eta(0);
    for (int i = 0; i < n_obs; i++) {
      mu(i) = theta(2 * i);
      sigma(i) = stan::math::exp(
          multiply(0.5, theta(2 * i + 1)));  // * sigma_global;
    }
    // return stan::math::normal_lpdf(y, mu, sigma);
    return stan::math::normal_lpdf(y, mu, multiply(sigma_global, sigma));
  }
};

TEST_F(laplace_motorcyle_gp_test, gp_motorcycle2) {
  using stan::math::laplace_marginal;
  using stan::math::laplace_marginal_tol;
  using stan::math::value_of;
  {
    double tolerance = 1e-12;
    constexpr int max_num_steps = 300;
    int hessian_block_size = 2;
    solver = 3;
    int do_line_search = 1;
    int max_steps_line_search = 10;
    double target = laplace_marginal_tol<false>(
        normal_likelihood2{}, std::forward_as_tuple(y, delta_int, eta),
        covariance_motorcycle_functor{},
        std::forward_as_tuple(x, phi_dbl(0), phi_dbl(1), phi_dbl(2), phi_dbl(3),
                              n_obs),
        theta0, tolerance, max_num_steps, hessian_block_size, solver,
        max_steps_line_search, nullptr);
  }
  // TODO(Charles): benchmark this result against GPStuff.
  constexpr double tolerance = 1e-8;
  constexpr int max_num_steps = 100;
  stan::test::ad_tolerances tols;
  tols.gradient_grad_ = 1e-3;

  stan::math::test::run_solver_grid(
      [&](int solver_num, int hessian_block_size, int max_steps_line_search,
          auto&& theta_0) {
        auto f = [&](auto&& eta_v, auto&& phi_0, auto&& phi) {
          return laplace_marginal_tol<false>(
              normal_likelihood2{}, std::forward_as_tuple(y, delta_int, eta_v),
              covariance_motorcycle_functor{},
              std::forward_as_tuple(x, phi_0, phi(1), phi(2), phi(3), n_obs),
              theta_0, tolerance, max_num_steps, hessian_block_size, solver_num,
              max_steps_line_search, nullptr);
        };
        stan::test::expect_ad<true>(tols, f, eta_dbl, phi_dbl(0), phi_dbl);
      },
      theta0);
}
