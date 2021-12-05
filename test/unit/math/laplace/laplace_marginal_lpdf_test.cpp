#include <stan/math.hpp>
#include <stan/math/laplace/laplace.hpp>
#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/rev/fun/util.hpp>
#include <stan/math/prim/fun/lgamma.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

struct poisson_log_likelihood2 {
  template <typename T_theta, typename T_eta>
  stan::return_type_t<T_theta, T_eta> operator() (
    const Eigen::Matrix<T_theta, -1, 1>& theta,
    const Eigen::Matrix<T_eta, -1, 1>& eta,
    const Eigen::VectorXd& y,
    const std::vector<int>& delta_int,
    std::ostream* pstream) const {
    return stan::math::poisson_log_lpmf(delta_int, theta);
  }
};

TEST(laplace_marginal_lpdf, poisson_log_phi_dim_2) {
  using stan::math::laplace_marginal_lpmf;
  using stan::math::to_vector;
  using stan::math::value_of;
  using stan::math::var;

  int dim_phi = 2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi(dim_phi);
  phi << 1.6, 0.45;

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

  Eigen::Matrix<var, Eigen::Dynamic, 1> eta_dummy;
  Eigen::VectorXd y_dummy;
  std::vector<double> delta_dummy;
  std::vector<int> delta_int_dummy;

  std::vector<int> n_samples = {1, 1};
  std::vector<int> sums = {1, 0};

  stan::math::test::squared_kernel_functor K;
  var target
    = laplace_marginal_lpmf<FALSE>(sums, poisson_log_likelihood2(),
                                   eta_dummy, y_dummy,
                                   K, phi, x, delta_dummy,
                                   delta_int_dummy, theta_0);

  // TODO: benchmark target against gpstuff.
  // Expected: -2.53056
  double tol = 1e-4;
  EXPECT_NEAR(-2.53056, value_of(target), tol);

  // Test with optional arguments
  double tolerance = 1e-6;
  int max_num_steps = 100;
  int hessian_block_size = 0;
  int solver = 1;
  int do_line_search = 1;
  int max_steps_line_search = 10;

  target = laplace_marginal_lpmf<FALSE>(sums, poisson_log_likelihood2(),
                                        eta_dummy, y_dummy,
                                        K, phi, x, delta_dummy,
                                        delta_int_dummy, theta_0,
                                        tolerance, max_num_steps,
                                        hessian_block_size, solver,
                                        do_line_search, max_steps_line_search);
  EXPECT_NEAR(-2.53056, value_of(target), tol);

  std::vector<double> g;
  std::vector<stan::math::var> parm_vec{phi(0), phi(1)};
  target.grad(parm_vec, g);

  // // finite diff test
  double diff = 1e-7;
  Eigen::VectorXd phi_dbl = value_of(phi);
  Eigen::VectorXd phi_1l = phi_dbl, phi_1u = phi_dbl, phi_2l = phi_dbl,
                  phi_2u = phi_dbl;
  phi_1l(0) -= diff;
  phi_1u(0) += diff;
  phi_2l(1) -= diff;
  phi_2u(1) += diff;

  Eigen::VectorXd eta_dummy_dbl = value_of(eta_dummy);

  double target_1u = laplace_marginal_lpmf<FALSE>(
    sums, poisson_log_likelihood2(), eta_dummy_dbl, y_dummy, K, phi_1u, x,
    delta_dummy, delta_int_dummy, theta_0),
  target_1l = laplace_marginal_lpmf<FALSE>(
    sums, poisson_log_likelihood2(), eta_dummy_dbl, y_dummy, K, phi_1l, x,
    delta_dummy, delta_int_dummy, theta_0),
  target_2u = laplace_marginal_lpmf<FALSE>(
    sums, poisson_log_likelihood2(), eta_dummy_dbl, y_dummy, K, phi_2u, x,
    delta_dummy, delta_int_dummy, theta_0),
  target_2l = laplace_marginal_lpmf<FALSE>(
    sums, poisson_log_likelihood2(), eta_dummy_dbl, y_dummy, K, phi_2l, x,
    delta_dummy, delta_int_dummy, theta_0);

  std::vector<double> g_finite(dim_phi);
  g_finite[0] = (target_1u - target_1l) / (2 * diff);
  g_finite[1] = (target_2u - target_2l) / (2 * diff);

  tol = 1.1e-4;
  EXPECT_NEAR(g_finite[0], g[0], tol);
  EXPECT_NEAR(g_finite[1], g[1], tol);
}

struct poisson_log_exposure_likelihood {
  template <typename T_theta, typename T_eta>
  stan::return_type_t<T_theta, T_eta> operator() (
    const Eigen::Matrix<T_theta, -1, 1>& theta,
    const Eigen::Matrix<T_eta, -1, 1>& eta,
    const Eigen::VectorXd& ye,
    const std::vector<int>& delta_int,
    std::ostream* pstream) const {
      return stan::math::poisson_log_lpmf(delta_int, theta
                                                       + stan::math::log(ye));
  }
};

TEST_F(laplace_disease_map_test, laplace_marginal_lpmf) {
  using stan::math::laplace_marginal_lpmf;
  using stan::math::laplace_marginal_poisson_log_lpmf;
  using stan::math::var;
  using stan::math::value_of;

  Eigen::Matrix<var, Eigen::Dynamic, 1> eta_dummy;
  // Eigen::VectorXd y_dummy;
  std::vector<double> delta_dummy;
  std::vector<int> delta_int_dummy;
  stan::math::test::sqr_exp_kernel_functor K;

  double tolerance = 1e-6;
  int max_num_steps = 100;
  int hessian_block_size = 0;
  int solver = 1;
  int do_line_search = 0;
  int max_steps_line_search = 10;

  var marginal_density
    = laplace_marginal_lpmf<FALSE>(y, poisson_log_exposure_likelihood(),
                                   eta_dummy, ye,
                                   K, phi, x, delta, delta_int, theta_0);

  double tol = 6e-4;
  // Benchmark from GPStuff.
  EXPECT_NEAR(-2866.88, value_of(marginal_density), tol);

  std::vector<double> g;
  std::vector<var> parm_vec{phi(0), phi(1)};
  marginal_density.grad(parm_vec, g);

  // finite diff
  Eigen::VectorXd phi_dbl = value_of(phi);
  Eigen::VectorXd phi_u0 = phi_dbl, phi_u1 = phi_dbl, phi_l0 = phi_dbl,
                  phi_l1 = phi_dbl;
  double eps = 1e-7;

  phi_u0(0) += eps;
  phi_u1(1) += eps;
  phi_l0(0) -= eps;
  phi_l1(1) -= eps;

  Eigen::VectorXd eta_dummy_dbl = value_of(eta_dummy);

  double target_u0 = laplace_marginal_lpmf<FALSE>(
      y, poisson_log_exposure_likelihood(), eta_dummy_dbl, ye, K, phi_u0, x,
      delta, delta_int, theta_0),
    target_u1 = laplace_marginal_lpmf<FALSE>(
        y, poisson_log_exposure_likelihood(), eta_dummy_dbl, ye, K, phi_u1, x,
        delta, delta_int, theta_0),
    target_l0 = laplace_marginal_lpmf<FALSE>(
            y, poisson_log_exposure_likelihood(), eta_dummy_dbl, ye, K, phi_l0, x,
          delta, delta_int, theta_0),
    target_l1 = laplace_marginal_lpmf<FALSE>(
              y, poisson_log_exposure_likelihood(), eta_dummy_dbl, ye, K, phi_l1, x,
              delta, delta_int, theta_0);

  EXPECT_NEAR((target_u0 - target_l0) / (2 * eps), g[0], 3e-3);
  EXPECT_NEAR((target_u1 - target_l1) / (2 * eps), g[1], 2e-4);
}

struct bernoulli_logit_likelihood {
  template <typename T_theta, typename T_eta>
  stan::return_type_t<T_theta, T_eta> operator() (
    const Eigen::Matrix<T_theta, -1, 1>& theta,
    const Eigen::Matrix<T_eta, -1, 1>& eta,
    const Eigen::VectorXd& ye,
    const std::vector<int>& delta_int,
    std::ostream* pstream) const {
      return stan::math::bernoulli_logit_lpmf(delta_int, theta);
  }
};

TEST(laplace_marginal_lpdf, bernoulli_logit_phi_dim500) {
  using stan::math::laplace_marginal_lpmf;
  using stan::math::to_vector;
  using stan::math::var;

  int dim_theta = 500;
  int n_observations = 500;
  std::string data_directory = "test/unit/math/laplace/aki_synth_data/";
  std::vector<double> x1(dim_theta), x2(dim_theta);
  std::vector<int> y(n_observations);
  stan::math::test::read_in_data(dim_theta, n_observations, data_directory, x1,
                                 x2, y);

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
  std::vector<int> delta_int;
  int dim_phi = 2;
  Eigen::Matrix<var, Eigen::Dynamic, 1> phi(dim_phi);
  phi << 1.6, 1;
  Eigen::Matrix<var, Eigen::Dynamic, 1> eta_dummy;

  stan::math::test::sqr_exp_kernel_functor K;
  bernoulli_logit_likelihood L;
  var target
    = laplace_marginal_lpmf<FALSE>(y, L, eta_dummy, delta_L,
                                   K, phi, x, delta, delta_int,
                                   theta_0);

  double tol = 8e-5;
  // Benchmark against gpstuff.
  EXPECT_NEAR(-195.368, value_of(target), tol);

  std::vector<double> g;
  std::vector<stan::math::var> parm_vec{phi(0), phi(1)};
  target.grad(parm_vec, g);

  // finite diff benchmark
  double diff = 1e-7;
  Eigen::VectorXd phi_dbl = value_of(phi);
  Eigen::VectorXd phi_1l = phi_dbl, phi_1u = phi_dbl, phi_2l = phi_dbl,
                  phi_2u = phi_dbl;
  phi_1l(0) -= diff;
  phi_1u(0) += diff;
  phi_2l(1) -= diff;
  phi_2u(1) += diff;
  Eigen::VectorXd eta_dummy_dbl;

  double target_1u = laplace_marginal_lpmf<FALSE>(y, L, eta_dummy_dbl, delta_L,
                                 K, phi_1u, x, delta, delta_int,
                                 theta_0),
         target_1l = laplace_marginal_lpmf<FALSE>(y, L, eta_dummy_dbl, delta_L,
                                        K, phi_1l, x, delta, delta_int,
                                        theta_0),
         target_2u = laplace_marginal_lpmf<FALSE>(y, L, eta_dummy_dbl, delta_L,
                                        K, phi_2u, x, delta, delta_int,
                                        theta_0),
         target_2l = laplace_marginal_lpmf<FALSE>(y, L, eta_dummy_dbl, delta_L,
                                        K, phi_2l, x, delta, delta_int,
                                        theta_0);

  std::vector<double> g_finite(dim_phi);
  g_finite[0] = (target_1u - target_1l) / (2 * diff);
  g_finite[1] = (target_2u - target_2l) / (2 * diff);

  tol = 1e-3;
  EXPECT_NEAR(g_finite[0], g[0], tol);
  EXPECT_NEAR(g_finite[1], g[1], tol);
}

struct covariance_motorcycle_functor {
  template <typename T1, typename T2>
  Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> operator()(
      const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi, const T2& x,
      const std::vector<double>& delta, const std::vector<int>& delta_int,
      std::ostream* msgs = nullptr) const {
    using Eigen::Matrix;
    using stan::math::gp_exp_quad_cov;

    T1 length_scale_f = phi(0);
    T1 length_scale_g = phi(1);
    T1 sigma_f = phi(2);
    T1 sigma_g = phi(3);
    int n_obs = delta_int[0];

    double jitter = 1e-6;
    Matrix<T1, -1, -1> kernel_f = gp_exp_quad_cov(x, sigma_f, length_scale_f);
    Matrix<T1, -1, -1> kernel_g = gp_exp_quad_cov(x, sigma_g, length_scale_g);

    Matrix<T1, -1, -1> kernel_all = Eigen::MatrixXd::Zero(2 * n_obs, 2 * n_obs);
    for (int i = 0; i < n_obs; i++) {
      for (int j = 0; j <= i; j++) {
        kernel_all(2 * i, 2 * j) = kernel_f(i, j);
        kernel_all(2 * i + 1, 2 * j + 1) = kernel_g(i, j);
        if (i != j) {
          kernel_all(2 * j, 2 * i) = kernel_all(2 * i, 2 * j);
          kernel_all(2 * j + 1, 2 * i + 1) = kernel_all(2 * i + 1, 2 * j + 1);
        }
      }
    }

    for (int i = 0; i < 2 * n_obs; i++)
      kernel_all(i, i) += jitter;

    return kernel_all;
  }
};

struct normal_likelihood {
  template <typename T_theta, typename T_eta>
  stan::return_type_t<T_theta, T_eta> operator()(
      const Eigen::Matrix<T_theta, -1, 1>& theta,
      const Eigen::Matrix<T_eta, -1, 1>& eta, const Eigen::VectorXd& y,
      const std::vector<int>& delta_int, std::ostream* pstream) const {
    int n_obs = delta_int[0];
    Eigen::Matrix<T_theta, -1, 1> mu(n_obs);
    Eigen::Matrix<T_theta, -1, 1> sigma(n_obs);
    for (int i = 0; i < n_obs; i++) {
      mu(i) = theta(2 * i);
      sigma(i) = exp(0.5 * theta(2 * i + 1));
    }

    return stan::math::normal_lpdf(y, mu, sigma);
  }
};

struct normal_likelihood2 {
  template <typename T_theta, typename T_eta>
  stan::return_type_t<T_theta, T_eta> operator()(
      const Eigen::Matrix<T_theta, -1, 1>& theta,
      const Eigen::Matrix<T_eta, -1, 1>& eta, const Eigen::VectorXd& y,
      const std::vector<int>& delta_int, std::ostream* pstream) const {
    using stan::math::multiply;
    int n_obs = delta_int[0];
    Eigen::Matrix<T_theta, -1, 1> mu(n_obs);
    Eigen::Matrix<T_theta, -1, 1> sigma(n_obs);
    T_eta sigma_global = eta(0);
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

    if (false) {
      n_obs = 6;
      Eigen::VectorXd x_vec(n_obs);
      x_vec << 2.4, 2.6, 3.2, 3.6, 4.0, 6.2;
      x.resize(n_obs);
      for (int i = 0; i < n_obs; i++)
        x[i] = x_vec(i);
      y.resize(n_obs);
      y << 0.0, -1.3, -2.7, 0.0, -2.7, -2.7;
    }

    if (true) {
      n_obs = 133;
      stan::math::test::read_data(
          n_obs, "test/unit/math/laplace/motorcycle_gp/", x, y);
    }

    length_scale_f = 0.3;
    length_scale_g = 0.5;
    sigma_f = 0.25;
    sigma_g = 0.25;

    dim_phi = 4;
    phi.resize(dim_phi);
    phi << length_scale_f, length_scale_g, sigma_f, sigma_g;

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
    phi_dbl = value_of(phi);
    eta_dbl = value_of(eta);
  }

  int n_obs;
  int dim_phi;
  std::vector<double> x;
  Eigen::VectorXd y;

  stan::math::var length_scale_f;
  stan::math::var length_scale_g;
  stan::math::var sigma_f;
  stan::math::var sigma_g;
  Eigen::Matrix<stan::math::var, -1, 1> phi;
  std::vector<int> delta_int;
  std::vector<double> delta_dummy;
  Eigen::VectorXd theta0;
  Eigen::VectorXd eta_dbl;
  Eigen::Matrix<stan::math::var, -1, 1> eta;
  int solver;
  double eps;
  Eigen::VectorXd phi_dbl;
};

TEST_F(laplace_motorcyle_gp_test, gp_motorcycle) {
  using stan::math::laplace_marginal_lpdf;
  using stan::math::value_of;
  using stan::math::var;

  covariance_motorcycle_functor K_f;
  normal_likelihood L;

  double tolerance = 1e-08;
  int max_num_steps = 100;
  int hessian_block_size = 2;
  solver = 3;
  int do_line_search = 1;
  int max_steps_line_search = 10;

  covariance_motorcycle_functor K;
  var target
    = laplace_marginal_lpdf<FALSE>(y, L, eta, delta_int,
                            K, phi, x, delta_dummy, delta_int, theta0,
                            tolerance, max_num_steps, hessian_block_size,
                            solver, do_line_search, max_steps_line_search);

  // TODO: benchmark this result against GPStuff.

  std::vector<double> g;
  std::vector<stan::math::var> parm_vec{phi(0), phi(1), phi(2), phi(3), eta(0)};
  target.grad(parm_vec, g);

  // finite diff benchmark
  double g_finite;
  for (int i = 0; i < dim_phi; i++) {
    Eigen::VectorXd phi_u = phi_dbl, phi_l = phi_dbl;
    phi_u(i) += eps;
    phi_l(i) -= eps;

    double target_u
      = laplace_marginal_lpdf<FALSE>(y, L, eta_dbl, delta_int,
                              K, phi_u, x, delta_dummy, delta_int, theta0,
                              tolerance, max_num_steps, hessian_block_size,
                              solver, do_line_search, max_steps_line_search);

    double target_l
      = laplace_marginal_lpdf<FALSE>(y, L, eta_dbl, delta_int,
                            K, phi_l, x, delta_dummy, delta_int, theta0,
                            tolerance, max_num_steps, hessian_block_size,
                            solver, do_line_search, max_steps_line_search);

    g_finite = (target_u - target_l) / (2 * eps);

    double tol = 1.2e-5;
    EXPECT_NEAR(g_finite, g[i], tol);
  }

  Eigen::VectorXd eta_u = eta_dbl, eta_l = eta_dbl;
  eta_u(0) += eps;
  eta_l(0) -= eps;

  double target_u
    = laplace_marginal_lpdf<FALSE>(y, L, eta_u, delta_int,
                            K, phi_dbl, x, delta_dummy, delta_int, theta0,
                            tolerance, max_num_steps, hessian_block_size,
                            solver, do_line_search, max_steps_line_search);

  double target_l
    = laplace_marginal_lpdf<FALSE>(y, L, eta_l, delta_int,
                            K, phi_dbl, x, delta_dummy, delta_int, theta0,
                            tolerance, max_num_steps, hessian_block_size,
                            solver, do_line_search, max_steps_line_search);

  g_finite = (target_u - target_l) / (2 * eps);
  double tol = 1e-7;
  EXPECT_NEAR(g_finite, g[dim_phi + 1], tol);
}
