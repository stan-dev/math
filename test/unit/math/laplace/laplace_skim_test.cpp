#include <stan/math.hpp>
#include <stan/math/laplace/laplace_likelihood.hpp>
#include <stan/math/laplace/laplace_marginal_bernoulli.hpp>

#include <test/unit/math/laplace/lgp_utility.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>

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

    // CHECK -- how does Stan create such an object?
    // Eigen::Matrix<T, -1, 1> one_matrix(N, N);
    // for (int i = 0; i < N; i++)
    //   for (int j = 0; j < N; j++) one_matrix(i, j) = 1;
    
    // std::cout << X.rows() << " " << X.cols() << std::endl
    //           << lambda_tilde.size() << std::endl
    //           << lambda_tilde.transpose() << std::endl; 
              // << diag_post_multiply(X, lambda_tilde).rows() << std::endl;

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

TEST(laplace, skm) {
  using stan::math::diff_logistic_log;
  using stan::math::var;
  using stan::math::square;
  using stan::math::elt_divide;
  using stan::math::add;
  using Eigen::MatrixXd;
  using Eigen::VectorXd;

  typedef Eigen::Matrix<var, -1, 1> Vector_v;
  typedef Eigen::Matrix<var, -1, -1> Matrix_v;

  // DATA AND TRANSFORMED DATA BLOCK
  int N = 100;
  int M = 2;  // options: 2, 50, 100, 150, 200

  std::string data_directory = "test/unit/math/laplace/skim_data/" +
    std::to_string(M) + "_" + std::to_string(N) + "/";
  MatrixXd X(N, M);
  std::vector<int> y(N);
  VectorXd lambda(M);

  read_in_data(M, N, data_directory, X, y, lambda);

  // std::cout << X << std::endl;
  // std::cout << lambda.transpose() << std::endl;
  // for (int i = 0; i < N; i++) std::cout << y[i] << " ";
  // std::cout << std::endl;

  double alpha_base = 0, psi = 1, m0 = 1,  // options: m0 = 2
    slab_scale = 3,
    slab_scale2 = slab_scale * slab_scale,
    slab_df = 25,
    half_slab_df = 0.5 * slab_df;

  VectorXd mu = VectorXd::Zero(N);
  std::vector<double> delta(1);
  delta[0] = psi;
  std::vector<int> delta_int(2);
  delta_int[0] = N;
  delta_int[1] = M;

  std::vector<int> n_samples(N, 1);
  VectorXd theta_0 = VectorXd::Zero(N);

  std::vector<VectorXd> x_tot(2 * N);
  {
    MatrixXd X2 = square(X);
    for (int n = 0; n < N; n++) x_tot[n] = X.block(n, 0, 1, M).transpose();
    for (int n = 0; n < N; n++) x_tot[N + n] = X2.block(n, 0, 1, M).transpose();
  }
  
  // PARAMETERS BLOCK
  // lambda term is defined above
  var c2_tilde = 1.112843,
    tau_tilde = 7.615908,
    sigma = 1.708423,
    eta_base = 0.9910583;
  
  // TRANSFORMED PARAMETERS BLOCK
  var phi = (m0 / (M - m0)) * (sigma / sqrt(N)) * tau_tilde,
    c2 = slab_scale2 * c2_tilde,
    eta = square(phi) / c2 * eta_base,
    alpha = square(phi) / c2 * alpha_base;
  
  Vector_v lambda_tilde = 
    c2 * elt_divide(square(lambda),
                    add(c2, multiply(square(phi), square(lambda))));

  Vector_v parm(M + 4);
  parm.head(M) = lambda_tilde;
  parm(M) = eta;
  parm(M + 1) = alpha;
  parm(M + 2) = phi;
  parm(M + 3) = sigma;

  K_functor K;
  // for (int i = 0; i < parm.size(); i++) std::cout << parm(i) << " ";
  // std::cout << std::endl;
  // std::cout << "x_tot" << std::endl;
  // for (size_t i = 0; i < x_tot.size(); i++) std::cout << x_tot[i].transpose() << std::endl;
  // std::cout << std::endl << std::endl;
  // for (size_t i = 0; i < delta.size(); i++) std::cout << delta[i] << std::endl;
  // for (size_t i = 0; i < delta_int.size(); i++) std::cout << delta_int[i] << std::endl;
    
  // std::cout << K(parm, x_tot, delta, delta_int, 0) << std::endl;

  auto start = std::chrono::system_clock::now();
  
  var marginal_density = laplace_marginal_bernoulli(y, n_samples, K_functor(),
                                      parm, x_tot, delta, delta_int, theta_0);
  
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;
  
  VEC g;
  AVEC parm_vec(M);
  for (int m = 0; m < M; m++) parm_vec[m] = parm(m);
  marginal_density.grad(parm_vec, g);

  std::cout << "LAPLACE MARGINAL AND VARI CLASS"
            << "density: " << marginal_density << std::endl
            << "autodiff grad: ";
  // for (size_t i = 0; i < parm.size(); i++) std::cout << g[i] << " ";
  std::cout << std::endl
            << "total time: " << elapsed_time.count() << std::endl
            << std::endl;
}
