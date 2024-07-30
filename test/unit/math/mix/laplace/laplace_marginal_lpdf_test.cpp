#include <test/unit/math/test_ad.hpp>
#include <stan/math.hpp>
#include <stan/math/mix.hpp>
#include <test/unit/math/mix/laplace/laplace_utility.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <test/unit/math/mix/laplace/aki_synth_data/x1.hpp>
#include <test/unit/math/mix/laplace/motorcycle_gp/x_vec.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

struct poisson_log_likelihood2 {
  template <typename Theta, typename Eta>
  auto operator()(const Theta& theta, const Eta& /* eta */,
                  const std::vector<int>& delta_int,
                  std::ostream* pstream) const {
    return stan::math::poisson_log_lpmf(delta_int, theta);
  }
};

TEST(laplace_marginal_lpdf, poisson_log_phi_dim_2) {
  using stan::math::laplace_marginal_lpmf;
  using stan::math::laplace_marginal_tol_lpmf;
  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;

  int dim_phi = 2;
  Eigen::Matrix<double, Eigen::Dynamic, 1> phi_dbl(dim_phi);
  phi_dbl << 1.6, 0.45;

  int dim_theta = 2;
  Eigen::VectorXd theta_0(dim_theta);
  theta_0 << 0, 0;

  int dim_x = 2;
  std::vector<Eigen::VectorXd> x(dim_theta);
  Eigen::VectorXd x_0(2);
  x_0 << 0.05100797, 0.16086164;
  Eigen::VectorXd x_1(2);
  x_1 << -0.59823393, 0.98701425;
  x[0] = x_0;
  x[1] = x_1;

  Eigen::Matrix<double, Eigen::Dynamic, 1> eta_dummy;
  Eigen::VectorXd y_dummy;
  std::vector<double> delta_dummy;
  std::vector<int> delta_int_dummy;

  std::vector<int> n_samples = {1, 1};
  std::vector<int> sums = {1, 0};

  stan::math::test::squared_kernel_functor K;
  double target = laplace_marginal_lpmf<false>(poisson_log_likelihood2(),
    std::forward_as_tuple(sums), theta_0, K, nullptr,
    x, phi_dbl(0), phi_dbl(1));

  // TODO: benchmark target against gpstuff.
  // Expected: -2.53056
  double tol = 1e-4;
  EXPECT_NEAR(-2.53056, value_of(target), tol);

  // Test with optional arguments
  {
  double tolerance = 1e-6;
  int max_num_steps = 100;
  int hessian_block_size = 1;
  int solver = 1;
  int do_line_search = 1;
  int max_steps_line_search = 10;

  target = laplace_marginal_tol_lpmf<false>(
      poisson_log_likelihood2(), std::forward_as_tuple(sums), theta_0,
      K, tolerance,
      max_num_steps, hessian_block_size, solver, max_steps_line_search, nullptr, x, phi_dbl(0), phi_dbl(1));
  EXPECT_NEAR(-2.53056, value_of(target), tol);
  }

  double tolerance = 1e-6;
  int max_num_steps = 100;
  stan::test::ad_tolerances ad_tol;
  ad_tol.gradient_val_ = 4e-4;
  ad_tol.gradient_grad_ = 1.1e-3;
  //FIXME(Steve): hessian_block_size of 3 fails approx test
  for (int max_steps_line_search = 0; max_steps_line_search < 4; ++max_steps_line_search) {
    for (int hessian_block_size = 1; hessian_block_size < 3; hessian_block_size++) {
      for (int solver_num = 1; solver_num < 4; solver_num++) {
        auto f = [&](auto&& alpha, auto&& rho) {
          return laplace_marginal_tol_lpmf<false>(
            poisson_log_likelihood2(), std::forward_as_tuple(sums), theta_0,
            K, tolerance, max_num_steps, hessian_block_size, solver_num,
            max_steps_line_search, nullptr, x, alpha, rho);
        };
        stan::test::expect_ad<true>(ad_tol, f, phi_dbl[0], phi_dbl[1]);
        }
    }
  }
}

struct poisson_log_exposure_likelihood {
  template <typename Theta, typename Eta>
  auto operator()(const Theta& theta, const Eta& /* eta */, const Eigen::VectorXd& ye,
                  const std::vector<int>& delta_int,
                  std::ostream* pstream) const {
    return stan::math::poisson_log_lpmf(delta_int, theta + stan::math::log(ye));
  }
};

TEST_F(laplace_disease_map_test, laplace_marginal_lpmf) {
  using stan::math::laplace_marginal_lpmf;
  using stan::math::laplace_marginal_tol_lpmf;
  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::value_of;
  using stan::math::var;

  Eigen::Matrix<double, -1, 1> eta_dummy;
  // Eigen::VectorXd y_dummy;
  std::vector<double> delta_dummy;
  std::vector<int> delta_int_dummy;
  stan::math::test::sqr_exp_kernel_functor K;

  double marginal_density = laplace_marginal_lpmf<false>(
      poisson_log_exposure_likelihood(), std::forward_as_tuple(ye, y), theta_0, K, nullptr,
      x, phi_dbl(0), phi_dbl(1));

  double tol = 6e-4;
  // Benchmark from GPStuff.
  EXPECT_NEAR(-2866.88, value_of(marginal_density), tol);
  double tolerance = 1e-6;
  int max_num_steps = 100;
  stan::test::ad_tolerances ad_tol;
  ad_tol.gradient_val_ = 8e-4;
  ad_tol.gradient_grad_ = 1.1e-3;
  //FIXME(Steve): hessian_block_size of 3 fails approx test
  for (int max_steps_line_search = 0; max_steps_line_search < 4; ++max_steps_line_search) {
    for (int hessian_block_size = 1; hessian_block_size < 3; hessian_block_size++) {
      for (int solver_num = 1; solver_num < 4; solver_num++) {
        auto f = [&](auto&& alpha, auto&& rho) {
          return laplace_marginal_tol_lpmf<false>(
            poisson_log_exposure_likelihood(), std::forward_as_tuple(ye, y), theta_0,
            K, tolerance, max_num_steps, hessian_block_size, solver_num,
            max_steps_line_search, nullptr, x, alpha, rho);
        };
        stan::test::expect_ad<true>(ad_tol, f, phi_dbl[0], phi_dbl[1]);
        }
    }
  }

}

struct bernoulli_logit_likelihood {
  template <typename Theta, typename Eta>
  auto operator()(const Theta& theta, const Eta& eta,
                  const std::vector<int>& delta_int,
                  std::ostream* pstream) const {
    return stan::math::bernoulli_logit_lpmf(delta_int, theta);
  }
};

TEST(laplace_marginal_lpdf, bernoulli_logit_phi_dim500) {
  using stan::math::laplace_marginal_lpmf;
  using stan::math::laplace_marginal_tol_lpmf;
  using stan::math::to_vector;

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
  Eigen::Matrix<double, Eigen::Dynamic, 1> eta_dummy;

  stan::math::test::sqr_exp_kernel_functor K;
  double target = laplace_marginal_lpmf<false>(bernoulli_logit_likelihood{},
    std::forward_as_tuple(y), theta_0,
    K, nullptr, x, phi_dbl(0), phi_dbl(1));

  double tol = 8e-5;
  // Benchmark against gpstuff.
  EXPECT_NEAR(-195.368, target, tol);

  double tolerance = 1e-6;
  int max_num_steps = 100;
  stan::test::ad_tolerances ad_tol;
  ad_tol.gradient_val_ = 4e-4;
  ad_tol.gradient_grad_ = 1.1e-3;
  //FIXME(Steve): hessian_block_size of 3 fails approx test
  for (int max_steps_line_search = 0; max_steps_line_search < 4; ++max_steps_line_search) {
    for (int hessian_block_size = 1; hessian_block_size < 3; hessian_block_size++) {
      for (int solver_num = 1; solver_num < 4; solver_num++) {
        auto f = [&](auto&& alpha, auto&& rho) {
          return laplace_marginal_tol_lpmf<false>(bernoulli_logit_likelihood{},
            std::forward_as_tuple(y),
            theta_0, K, tolerance,
            max_num_steps, hessian_block_size, solver_num, max_steps_line_search, nullptr, x, alpha, rho);
        };
        stan::test::expect_ad<true>(ad_tol, f, phi_dbl[0], phi_dbl[1]);
        }
    }
  }
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

    double jitter = 1e-6;
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
  template <typename Theta, typename Eta>
  auto operator()(const Theta& theta, const Eta& eta, const Eigen::VectorXd& y,
                  const std::vector<int>& delta_int,
                  std::ostream* pstream) const {
    int n_obs = delta_int[0];
    Eigen::Matrix<stan::return_type_t<Theta>, -1, 1> mu(n_obs);
    Eigen::Matrix<stan::return_type_t<Theta>, -1, 1> sigma(n_obs);
    for (Eigen::Index i = 0; i < n_obs; i++) {
      mu(i) = theta(2 * i);
      sigma(i) = exp(0.5 * theta(2 * i + 1));
    }
    return stan::math::normal_lpdf(y, mu, sigma);
  }
};

struct normal_likelihood2 {
  template <typename Theta, typename Eta>
  auto operator()(const Theta& theta, const Eta& eta, const Eigen::VectorXd& y,
                  const std::vector<int>& delta_int,
                  std::ostream* pstream) const {
    using stan::math::multiply;
    int n_obs = delta_int[0];
    Eigen::Matrix<stan::return_type_t<Theta>, -1, 1> mu(n_obs);
    Eigen::Matrix<stan::return_type_t<Theta>, -1, 1> sigma(n_obs);
    stan::return_type_t<Eta> sigma_global = eta(0);
    for (int i = 0; i < n_obs; i++) {
      mu(i) = theta(2 * i);
      sigma(i) = exp(0.5 * theta(2 * i + 1));  // * sigma_global;
    }

    // return stan::math::normal_lpdf(y, mu, sigma);
    return stan::math::normal_lpdf(y, mu, multiply(sigma_global, sigma));
  }
};

class laplace_motorcyle_gp_test : public ::testing::Test {
 protected:
  void SetUp() override {
    using stan::math::gp_exp_quad_cov;
    using stan::math::value_of;



      n_obs = 133;
    x = stan::test::laplace::moto::x;
    y = stan::test::laplace::moto::y;

    length_scale_f = 0.3;
    length_scale_g = 0.5;
    sigma_f = 0.25;
    sigma_g = 0.25;

    dim_phi = 4;
    phi_dbl.resize(dim_phi);
    phi_dbl << length_scale_f, length_scale_g, sigma_f, sigma_g;

    eta.resize(1);
    eta(0) = 1;

    delta_int.resize(1);
    delta_int[0] = n_obs;

    theta0 = Eigen::VectorXd::Zero(2 * n_obs);
    // theta0 << -10, 0, -10, 0, -10, 0, -10,
    //           0, -10, 0, -10, 0;

    Eigen::MatrixXd K_plus_I
        = gp_exp_quad_cov(x, value_of(sigma_f), value_of(length_scale_f))
          + Eigen::MatrixXd::Identity(n_obs, n_obs);

    Eigen::VectorXd mu_hat = K_plus_I.colPivHouseholderQr().solve(y);

    // Remark: finds optimal point with or without informed initial guess.
    for (int i = 0; i < n_obs; i++) {
      theta0(2 * i) = mu_hat(i);  // 0
      theta0(2 * i + 1) = 0;
    }

    solver = 2;
    eps = 1e-7;
    eta_dbl = value_of(eta);
  }

  int n_obs;
  int dim_phi;
  std::vector<double> x;
  Eigen::VectorXd y;

  double length_scale_f;
  double length_scale_g;
  double sigma_f;
  double sigma_g;
  std::vector<int> delta_int;
  std::vector<double> delta_dummy;
  Eigen::VectorXd theta0;
  Eigen::VectorXd eta_dbl;
  Eigen::Matrix<double, -1, 1> eta;
  int solver;
  double eps;
  Eigen::VectorXd phi_dbl;
};

TEST_F(laplace_motorcyle_gp_test, gp_motorcycle) {
  using stan::math::laplace_marginal_lpdf;
  using stan::math::laplace_marginal_tol_lpdf;
  using stan::math::value_of;

  covariance_motorcycle_functor K_f;
  normal_likelihood L;
  {
  double tolerance = 1e-08;
  int max_num_steps = 100;
  int hessian_block_size = 2;
  solver = 3;
  int do_line_search = 1;
  int max_steps_line_search = 10;

  covariance_motorcycle_functor K;
  double target = laplace_marginal_tol_lpdf<false>(
      L, std::forward_as_tuple(y, delta_int), eta, theta0, K,
      tolerance, max_num_steps, hessian_block_size,
      solver, max_steps_line_search, nullptr, x, phi_dbl(0), phi_dbl(1),
      phi_dbl(2), phi_dbl(3), n_obs);
  }
  // TODO: benchmark this result against GPStuff.
  double tolerance = 1e-6;
  int max_num_steps = 100;
  stan::test::ad_tolerances ad_tol;
  ad_tol.gradient_val_ = 4e-4;
  ad_tol.gradient_grad_ = 1.1e-3;
  covariance_motorcycle_functor K;
  //FIXME(Steve): hessian_block_size of 3 fails approx test
  for (int max_steps_line_search = 0; max_steps_line_search < 4; ++max_steps_line_search) {
    for (int hessian_block_size = 1; hessian_block_size < 3; hessian_block_size++) {
      for (int solver_num = 1; solver_num < 4; solver_num++) {
        auto f = [&](auto&& phi) {
          return laplace_marginal_tol_lpdf<false>(
      L, std::forward_as_tuple(y, delta_int), eta, theta0, K,
      tolerance, max_num_steps, hessian_block_size,
      solver, max_steps_line_search, nullptr, x, phi(0), phi(1),
      phi(2), phi(3), n_obs);
        };
        stan::test::expect_ad<true>(ad_tol, f, phi_dbl);
        }
    }
  }

}
