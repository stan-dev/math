#include <stan/math.hpp>
#include <stan/math/laplace/laplace.hpp>
#include <stan/math/laplace/laplace_likelihood_general.hpp>
#include <stan/math/laplace/laplace_marginal_lpdf.hpp>

#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

// Example from
// https://avehtari.github.io/casestudies/Motorcycle/motorcycle.html

struct covariance_motorcycle_functor {
  template <typename T1, typename T2>
  Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic>
  operator() (const Eigen::Matrix<T1, Eigen::Dynamic, 1>& phi,
              const T2& x,
              const std::vector<double>& delta,
              const std::vector<int>& delta_int,
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

    Matrix<T1, -1, -1> kernel_all
     = Eigen::MatrixXd::Zero(2 * n_obs, 2 * n_obs);
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

    for (int i = 0; i < 2 * n_obs; i++) kernel_all(i, i) += jitter;

    return kernel_all;
  }
};

struct normal_likelihood {
  template<typename T_theta, typename T_eta>
  stan::return_type_t<T_theta, T_eta>
  operator()(const Eigen::Matrix<T_theta, -1, 1>& theta,
             const Eigen::Matrix<T_eta, -1, 1>& eta,
             const Eigen::VectorXd& y,
             const std::vector<int>& delta_int,
             std::ostream* pstream) const {
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

// include a global variance (passed through eta)
struct normal_likelihood2 {
  template<typename T_theta, typename T_eta>
  stan::return_type_t<T_theta, T_eta>
  operator()(const Eigen::Matrix<T_theta, -1, 1>& theta,
             const Eigen::Matrix<T_eta, -1, 1>& eta,
             const Eigen::VectorXd& y,
             const std::vector<int>& delta_int,
             std::ostream* pstream) const {
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

class laplace_motorcyle_gp_test : public::testing::Test {
protected:
  void SetUp() override {
    using stan::math::value_of;
    using stan::math::gp_exp_quad_cov;

    if (FALSE) {
      n_obs = 6;
      Eigen::VectorXd x_vec(n_obs);
      x_vec << 2.4, 2.6, 3.2, 3.6, 4.0, 6.2;
      x.resize(n_obs);
      for (int i = 0; i < n_obs; i++) x[i] = x_vec(i);
      y.resize(n_obs);
      y << 0.0, -1.3, -2.7,  0.0, -2.7, -2.7;
    }

    if (TRUE) {
      n_obs = 133;
      read_data(n_obs, "test/unit/math/laplace/motorcycle_gp/",
                   x, y);
      // std::cout << "x: ";
      // for (int i = 0; i < n_obs; i++) std::cout << x[i] << " ";
      // std::cout << " ..." << std::endl;
      // std::cout << "y: " << y.transpose().head(n_obs) << " ..." << std::endl;
    }

    // [0.335852,0.433641,0.335354,0.323559]
    length_scale_f = 0.335852;  // 0.3;
    length_scale_g = 0.433641;  // 0.5;
    sigma_f = 0.335354;  // 0.25;
    sigma_g = 0.323559;  // 0.25;

    phi.resize(4);
    phi << length_scale_f, length_scale_g, sigma_f, sigma_g;

    delta_int.resize(1);
    delta_int[0] = n_obs;

    theta0 = Eigen::VectorXd::Zero(2 * n_obs);
    // theta0 << -10, 0, -10, 0, -10, 0, -10,
    //           0, -10, 0, -10, 0;

    Eigen::MatrixXd
      K_plus_I = gp_exp_quad_cov(x, value_of(sigma_f), value_of(length_scale_f))
        + Eigen::MatrixXd::Identity(n_obs, n_obs);

    Eigen::VectorXd mu_hat
      = K_plus_I.colPivHouseholderQr().solve(y);

    // Remark: finds optimal point with or without informed initial guess.
    for (int i = 0; i < n_obs; i++) {
      theta0(2 * i) = mu_hat(i);  // 0
      theta0(2 * i + 1) = 0;
    }

    solver = 2;
  }

  int n_obs;
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
  Eigen::VectorXd eta_dummy_dbl;
  Eigen::Matrix<stan::math::var, -1, 1> eta_dummy;
  int solver;
};

TEST_F(laplace_motorcyle_gp_test, lk_autodiff) {
  using stan::math::var;
  using stan::math::value_of;
  using stan::math::laplace_marginal_density;
  using stan::math::diff_likelihood;

  covariance_motorcycle_functor K_f;
  Eigen::VectorXd phi_dbl_ = value_of(phi);
  Eigen::MatrixXd K_eval
    = K_f(phi_dbl_, x, delta_dummy, delta_int, 0);
  std::cout << "K_eval: " << K_eval.row(0).head(5) << std::endl;

  normal_likelihood f;
  diff_likelihood<normal_likelihood> diff_functor(f, y, delta_int);

  int hessian_block_size = 2;
  solver = 2;
  int do_line_search = 1;
  int max_steps_line_search = 100;
  double marginal_density_dbl
    = laplace_marginal_density(diff_functor,
                               covariance_motorcycle_functor(),
                               value_of(phi), eta_dummy_dbl,
                               x, delta_dummy, delta_int, theta0,
                               0, 1e-2, 20000, hessian_block_size,
                               solver, do_line_search,
                               max_steps_line_search);

  std::cout << "density: " << marginal_density_dbl << std::endl;

/*
  var marginal_density
    = laplace_marginal_density(diff_functor,
                               covariance_motorcycle_functor(),
                               phi, eta_dummy,
                               x, delta_dummy, delta_int, theta0,
                               0, 1e-8, 1000, hessian_block_size,
                               compute_W_root);

  VEC g;
  AVEC parm_vec = createAVEC(phi(0), phi(1), phi(2), phi(3));
  marginal_density.grad(parm_vec, g);
  std::cout << "grad: " << g[0] << " " << g[1] << " " << g[2] << " " << g[3]
  << std::endl;

  // FINITE DIFF benchmark
  // TODO: benchmark against all the inputs.
  double eps = 1e-7;
  Eigen::VectorXd phi_dbl = value_of(phi);
  Eigen::VectorXd phi_u0 = phi_dbl, phi_l0 = phi_dbl;
  phi_u0(0) += eps;
  phi_l0(0) -= eps;

  // TODO: test all the gradients in code (rather than doing it manually)
  double target_u0
    = laplace_marginal_density(diff_functor,
                               covariance_motorcycle_functor(),
                               phi_u0, eta_dummy_dbl,
                               x, delta_dummy, delta_int, theta0,
                               0, 1e-6, 100, hessian_block_size,
                               compute_W_root);


  double target_l0
    = laplace_marginal_density(diff_functor,
                               covariance_motorcycle_functor(),
                               phi_l0, eta_dummy_dbl,
                               x, delta_dummy, delta_int, theta0,
                               0, 1e-6, 100, hessian_block_size,
                               compute_W_root);

  std::cout << "g[0]: " << (target_u0 - target_l0) / (2 * eps) << std::endl;
  */
}

/*
TEST_F(laplace_motorcyle_gp_test, lk_autodiff_eta) {
  using stan::math::var;
  using stan::math::value_of;
  using stan::math::laplace_marginal_density;
  using stan::math::diff_likelihood;

  normal_likelihood2 f;
  diff_likelihood<normal_likelihood2> diff_functor(f, y, delta_int);
  Eigen::Matrix<var, -1, 1> eta(1);
  eta(0) = 1;
  int hessian_block_size = 2;
  double marginal_density_dbl
    = laplace_marginal_density(diff_functor,
                               covariance_motorcycle_functor(),
                               value_of(phi), value_of(eta),
                               x, delta_dummy, delta_int, theta0,
                               0, 1e-8, 100, hessian_block_size,
                               compute_W_root);

  std::cout << "density: " << marginal_density_dbl << std::endl;

  var marginal_density
    = laplace_marginal_density(diff_functor,
                               covariance_motorcycle_functor(),
                               phi, eta,
                               x, delta_dummy, delta_int, theta0,
                               0, 1e-8, 100, hessian_block_size,
                               compute_W_root);

  VEC g;
  AVEC parm_vec = createAVEC(phi(0), phi(1), phi(2), phi(3), eta(0));
  marginal_density.grad(parm_vec, g);
  std::cout << "grad: "
    << g[0] << " " << g[1] << " " << g[2] << " " << g[3] << " " << g[4]
    << std::endl;

  // finite diff benchmark
  // TODO: test all the gradients in code (rather than doing it manually)
  double eps = 1e-7;
  Eigen::VectorXd eta_dbl = value_of(eta);
  Eigen::VectorXd eta_u = eta_dbl, eta_l = eta_dbl;
  eta_u(0) += eps;
  eta_l(0) -= eps;

  double target_u
  = laplace_marginal_density(diff_functor,
                             covariance_motorcycle_functor(),
                             value_of(phi), eta_u,
                             x, delta_dummy, delta_int, theta0,
                             0, 1e-8, 100, hessian_block_size,
                             compute_W_root);

  double target_l
  = laplace_marginal_density(diff_functor,
                             covariance_motorcycle_functor(),
                             value_of(phi), eta_l,
                             x, delta_dummy, delta_int, theta0,
                             0, 1e-8, 100, hessian_block_size,
                             compute_W_root);

  std::cout << "gf[4]: " << (target_u - target_l) / (2 * eps) << std::endl;

} */
/*
TEST_F(laplace_motorcyle_gp_test, wrapper_function) {
  using stan::math::var;
  using stan::math::laplace_marginal_lpdf;

  // TODO: move this to the class test.
  Eigen::Matrix<var, -1, 1> eta(1);
  eta(0) = 1;
  int hessian_block_size = 2;

  var marginal_density
    = laplace_marginal_lpdf<FALSE>(y, normal_likelihood2(), eta, delta_int,
                            covariance_motorcycle_functor(), phi,
                            x, delta_dummy, delta_int, theta0,
                            1e-8, 100, hessian_block_size,
                            compute_W_root, 0);

  std::cout << "density: " << marginal_density << std::endl;

  VEC g;
  AVEC parm_vec = createAVEC(phi(0), phi(1), phi(2), phi(3), eta(0));
  marginal_density.grad(parm_vec, g);
  std::cout << "grad: "
  << g[0] << " " << g[1] << " " << g[2] << " " << g[3] << " " << g[4]
  << std::endl;
}  */
