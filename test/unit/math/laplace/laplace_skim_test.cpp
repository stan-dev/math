#include <stan/math.hpp>
#include <stan/math/laplace/laplace_likelihood_general.hpp>
#include <stan/math/laplace/laplace_likelihood_bernoulli_logit.hpp>
#include <stan/math/laplace/laplace_marginal_bernoulli_logit.hpp>

#include <test/unit/math/laplace/laplace_utility.hpp>
#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

struct K_functor {
  template <typename T>
  Eigen::Matrix<T, -1, -1>
  operator()(const Eigen::Matrix<T, -1, 1>& parm,
             const std::vector<Eigen::VectorXd>& x_tot,
             const std::vector<double>& delta,
             const std::vector<int>& delta_int,
             std::ostream* pstream) const {
    using stan::math::add;
    using stan::math::multiply;
    using stan::math::diag_post_multiply;
    using stan::math::square;
    using stan::math::transpose;

    int N = delta_int[0];
    int M = delta_int[1];

    Eigen::Matrix<T, -1, 1> lambda_tilde(M);
    for (int m = 0; m < M; m++) lambda_tilde[m] = parm[m];

    T eta = parm[M];
    T alpha = parm[M + 1];
    T phi = parm[M + 2];
    T sigma = parm[M + 3];
    double psi = delta[0];

    // Note -- when the object is declared as a scalar matrix,
    // the differentiation slows down.
    // Eigen::Matrix<T, -1, -1> X(N, M);
    // Eigen::Matrix<T, -1, -1> X2(N, M);
    Eigen::MatrixXd X(N, M);
    Eigen::MatrixXd X2(N, M);

    // CHECK -- does Stan really do a for loop here?
    // CHECK -- does Stan construct X and X2 as matrices of scalars
    //          and does this have an effect?
    // TO DO -- test when this is constructed using the block method.
    for (int n = 0; n < N; n++)
      for (int m = 0; m < M; m++) {
        X(n, m) = x_tot[n](m);
        X2(n, m) = x_tot[N + n](m);
      }

    Eigen::Matrix<T, -1, -1>
      K1 = multiply(diag_post_multiply(X, lambda_tilde), transpose(X));
    Eigen::Matrix<T, -1, -1>
      K2 = multiply(diag_post_multiply(X2, lambda_tilde), transpose(X2));

    Eigen::Matrix<T, -1, -1> K;
    K = square(eta) * square(add(K1, 1)) +
        (square(alpha) - 0.5 * square(eta)) * K2 +
        (square(phi) - square(eta)) * K1;
    K = add(0.5 + square(psi) - 0.5 * square(eta), K);

    // Add jitter to make linear algebra more numerically stable
    for (int n = 0; n < N; n++) K(n, n) += square(sigma) + 1e-7;
    return K;
  }
};

// Overload structure for case where x is passed as a matrix.
struct K_functor2 {
  template <typename T>
  Eigen::Matrix<T, -1, -1>
  operator()(const Eigen::Matrix<T, -1, 1>& parm,
             const Eigen::MatrixXd& x_tot,
             const std::vector<double>& delta,
             const std::vector<int>& delta_int,
             std::ostream* pstream) const {
    using stan::math::add;
    using stan::math::multiply;
    using stan::math::diag_post_multiply;
    using stan::math::square;
    using stan::math::transpose;

    int N = delta_int[0];
    int M = delta_int[1];

    Eigen::Matrix<T, -1, 1> lambda_tilde(M);
    for (int m = 0; m < M; m++) lambda_tilde[m] = parm[m];

    T eta = parm[M];
    T alpha = parm[M + 1];
    T phi = parm[M + 2];
    T sigma = parm[M + 3];
    double psi = delta[0];

    Eigen::MatrixXd X = x_tot.block(0, 0, N, M);
    Eigen::MatrixXd X2 = x_tot.block(N, 0, N, M);

    Eigen::Matrix<T, -1, -1>
      K1 = multiply(diag_post_multiply(X, lambda_tilde), transpose(X));
    Eigen::Matrix<T, -1, -1>
      K2 = multiply(diag_post_multiply(X2, lambda_tilde), transpose(X2));

    Eigen::Matrix<T, -1, -1> K;
    K = square(eta) * square(add(K1, 1)) +
        (square(alpha) - 0.5 * square(eta)) * K2 +
        (square(phi) - square(eta)) * K1;
    K = add(0.5 + square(psi) - 0.5 * square(eta), K);

    // Add jitter to make linear algebra more numerically stable
    for (int n = 0; n < N; n++) K(n, n) += square(sigma) + 1e-7;
    return K;
  }
};

class laplace_skim_test : public::testing::Test {
protected:
  void SetUp() override {
    using stan::math::square;
    using stan::math::var;
    using stan::math::square;
    using stan::math::elt_divide;
    using stan::math::add;

    N = 100;
    M = 200;  // options: 2, 50, 100, 150, 200
    // TODO: add to GitHub directory simulation for each configuration.
    // std::string data_directory = "test/unit/math/laplace/skim_data/" +
    //   std::to_string(M) + "_" + std::to_string(N) + "/";
    std::string data_directory = "test/unit/math/laplace/skim_data/";

    X.resize(N, M);
    y.resize(N);
    lambda.resize(M);

    read_in_data(M, N, data_directory, X, y, lambda);

    if (FALSE){
      std::cout << X << std::endl << "-----" << std::endl;
      std::cout << lambda.transpose() << std::endl << "------" << std::endl;
      std::cout << y[0] << " " << y[1] << " " << std::endl
        << "------" << std::endl;
    }

    alpha_base = 0;
    psi = 1;
    m0 = 1;
    slab_scale = 3;
    slab_scale2 = slab_scale * slab_scale;
    half_slab_df = 0.5 * slab_df;

    mu = Eigen::VectorXd::Zero(N);
    delta.resize(1);
    delta[0] = psi;
    delta_int.resize(2);
    delta_int[0] = N;
    delta_int[1] = M;

    std::vector<int> n_samples_(N, 1);
    n_samples = n_samples_;

    theta_0 = Eigen::VectorXd::Zero(N);

    X2 = square(X);
    x_tot_m.resize(2 * N, M);
    x_tot_m.block(0, 0, N, M) = X;
    x_tot_m.block(N, 0, N, M) = X2;

    // parameters block
    c2_tilde = 1.112843;
    tau_tilde = 7.615908;
    sigma = 1.708423;
    eta_base = 0.9910583;

    phi = (m0 / (M - m0)) * (sigma / sqrt(N)) * tau_tilde;
    c2 = slab_scale2 * c2_tilde;
    eta = square(phi) / c2 * eta_base;
    alpha = square(phi) / c2 * alpha_base;

    lambda_tilde = c2 * elt_divide(square(lambda),
                    add(c2, multiply(square(phi), square(lambda))));

    parm.resize(M + 4);
    parm.head(M) = lambda_tilde;
    parm(M) = eta;
    parm(M + 1) = alpha;
    parm(M + 2) = phi;
    parm(M + 3) = sigma;

    // std::cout << "parm: " << parm << std::endl;
  }

  int N;
  int M;
  Eigen::MatrixXd X;
  std::vector<int> y;
  Eigen::VectorXd lambda;
  double alpha_base, psi, m0, slab_scale, slab_scale2, slab_df, half_slab_df;
  Eigen::VectorXd mu;
  std::vector<double> delta;
  std::vector<int> delta_int;
  std::vector<int> n_samples;
  Eigen::VectorXd theta_0;
  Eigen::MatrixXd X2;
  Eigen::MatrixXd x_tot_m;

  stan::math::var c2_tilde, tau_tilde, sigma, eta_base,
    phi, c2, eta, alpha;
  Eigen::Matrix<stan::math::var, -1, 1> lambda_tilde;
  Eigen::Matrix<stan::math::var, -1, 1> parm;
};


TEST_F(laplace_skim_test, lk_analytical) {
  using stan::math::var;
  using stan::math::laplace_marginal_bernoulli_logit_lpmf;


  auto start = std::chrono::system_clock::now();

  var marginal_density
    = laplace_marginal_bernoulli_logit_lpmf(y, n_samples, K_functor2(),
                                            parm, x_tot_m, delta, delta_int,
                                            theta_0);

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;

  VEC g;
  AVEC parm_vec(M + 4);
  for (int m = 0; m < M + 4; m++) parm_vec[m] = parm(m);
  marginal_density.grad(parm_vec, g);

  // std::cout << parm << std::endl;

  // Expected density: - 10.9795
  std::cout << "LAPLACE MARGINAL AND VARI CLASS" << std::endl
            << "M: " << M << std::endl
            << "density: " << marginal_density << std::endl
            << "autodiff grad: ";
  for (size_t i = 0; i < 10; i++) std::cout << g[i] << " ";
  std::cout << std::endl
            << "total time: " << elapsed_time.count() << std::endl
            << std::endl;
}

struct bernoulli_logit_likelihood {
  template<typename T_theta, typename T_eta>
  stan::return_type_t<T_theta, T_eta>
  operator()(const Eigen::Matrix<T_theta, -1, 1>& theta,
             const Eigen::Matrix<T_eta, -1, 1>& eta,
             const Eigen::VectorXd& sums,       // sums
             const std::vector<int>& n_samples,  // n_samples
             std::ostream* pstream) const {
    using stan::math::to_vector;
    stan::math::diff_bernoulli_logit
      diff_functor(to_vector(n_samples), sums);

    return diff_functor.log_likelihood(theta, eta);
  }
};


TEST_F(laplace_skim_test, lk_autodiff) {
  using stan::math::var;
  using stan::math::laplace_marginal_density;
  using stan::math::diff_likelihood;
  using stan::math::to_vector;
  using stan::math::value_of;

  bernoulli_logit_likelihood f;
  diff_likelihood<bernoulli_logit_likelihood>
    diff_functor(f, to_vector(y), n_samples);

  auto start = std::chrono::system_clock::now();

  Eigen::Matrix<var, -1, 1> eta_dummy;
  var marginal_density
    = laplace_marginal_density(diff_functor, K_functor2(), parm, eta_dummy,
                               x_tot_m, delta, delta_int, theta_0);

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;

  VEC g;
  AVEC parm_vec(M + 4);
  for (int m = 0; m < M + 4; m++) parm_vec[m] = parm(m);
  marginal_density.grad(parm_vec, g);

  // Expected density: - 10.9795
  std::cout << "LAPLACE MARGINAL AND VARI CLASS" << std::endl
            << "M: " << M << std::endl
            << "density: " << marginal_density << std::endl
            << "autodiff grad: ";
  for (size_t i = 0; i < 10; i++) std::cout << g[i] << " ";
  std::cout << std::endl
            << "total time: " << elapsed_time.count() << std::endl
            << std::endl;
}
