#include <stan/math.hpp>
#include <stan/math/laplace/laplace_likelihood.hpp>
#include <stan/math/laplace/laplace_marginal.hpp>
#include <test/unit/math/rev/fun/util.hpp>

#include <gtest/gtest.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <vector>

TEST(laplace, likelihood_differentiation) {
  using stan::math::diff_neg_binomial_2_log;
  using stan::math::var;

  Eigen::VectorXd theta(2);
  theta << 1, 1;
  int n_theta = theta.size();
  Eigen::VectorXd eta(1);
  eta << 1.2;

  Eigen::VectorXd y(2);
  y << 0, 1;
  std::vector<int> y_index(2);
  y_index[0] = 0;
  y_index[1] = 1;

  // Eigen::Matrix<var, Eigen::Dynamic, 1> theta_v = theta;
  diff_neg_binomial_2_log diff_functor(y, y_index, n_theta);

  double log_density = diff_functor.log_likelihood(theta, eta);

  // benchmark against R
  EXPECT_FLOAT_EQ(-3.023328, log_density);

  Eigen::VectorXd gradient, hessian;
  diff_functor.diff(theta, eta, gradient, hessian);

  Eigen::VectorXd third_diff = diff_functor.third_diff(theta, eta);

  // Benchmark against finite diff
  double epsilon = 1e-6;
  Eigen::VectorXd theta_l0 = theta, theta_u0 = theta,
                  theta_l1 = theta, theta_u1 = theta;
  theta_u0(0) += epsilon;
  theta_l0(0) -= epsilon;
  theta_u1(1) += epsilon;
  theta_l1(1) -= epsilon;

  Eigen::VectorXd finite_gradient(2);
  finite_gradient(0) =
    (diff_functor.log_likelihood(theta_u0, eta)
      - diff_functor.log_likelihood(theta_l0, eta)) / (2 * epsilon);

  finite_gradient(1) =
    (diff_functor.log_likelihood(theta_u1, eta)
      - diff_functor.log_likelihood(theta_l1, eta)) / (2 * epsilon);

  Eigen::VectorXd gradient_l0, gradient_u0, gradient_l1, gradient_u1;
  Eigen::VectorXd hessian_l0, hessian_u0, hessian_l1, hessian_u1;
  Eigen::VectorXd hessian_dummy;
  diff_functor.diff(theta_l0, eta, gradient_l0, hessian_l0);
  diff_functor.diff(theta_u0, eta, gradient_u0, hessian_u0);
  diff_functor.diff(theta_l1, eta, gradient_l1, hessian_l1);
  diff_functor.diff(theta_u1, eta, gradient_u1, hessian_u1);

  Eigen::VectorXd finite_hessian(2);
  finite_hessian(0) = (gradient_u0 - gradient_l0)(0) / (2 * epsilon);
  finite_hessian(1) = (gradient_u1 - gradient_l1)(1) / (2 * epsilon);

  Eigen::VectorXd finite_third_diff(2);
  finite_third_diff(0) = (hessian_u0 - hessian_l0)(0) / (2 * epsilon);
  finite_third_diff(1) = (hessian_u1 - hessian_l1)(1) / (2 * epsilon);


  EXPECT_FLOAT_EQ(finite_gradient(0), gradient(0));
  EXPECT_FLOAT_EQ(finite_gradient(1), gradient(1));
  EXPECT_FLOAT_EQ(finite_hessian(0), hessian(0));
  EXPECT_FLOAT_EQ(finite_hessian(1), hessian(1));
  EXPECT_FLOAT_EQ(finite_third_diff(0), third_diff(0));
  EXPECT_FLOAT_EQ(finite_third_diff(1), third_diff(1));

  // derivatives wrt eta
  Eigen::VectorXd diff_eta = diff_functor.diff_eta(theta, eta);

  Eigen::VectorXd eta_l(1), eta_u(1);
  eta_l(0) = eta(0) - epsilon;
  eta_u(0) = eta(0) + epsilon;
  double finite_gradient_eta =
    (diff_functor.log_likelihood(theta, eta_u)
      - diff_functor.log_likelihood(theta, eta_l)) / (2 * epsilon);

  EXPECT_FLOAT_EQ(finite_gradient_eta,  diff_eta(0));

  Eigen::MatrixXd diff_theta_eta = diff_functor.diff_theta_eta(theta, eta);

  Eigen::VectorXd gradient_theta_l,
                  gradient_theta_u,
                  hessian_theta_u,
                  hessian_theta_l;

  diff_functor.diff(theta, eta_l, gradient_theta_l, hessian_theta_l);
  diff_functor.diff(theta, eta_u, gradient_theta_u, hessian_theta_u);
  Eigen::VectorXd finite_gradient_theta_eta
    = (gradient_theta_u - gradient_theta_l) / (2 * epsilon);

  EXPECT_FLOAT_EQ(finite_gradient_theta_eta(0), diff_theta_eta(0, 0));
  EXPECT_FLOAT_EQ(finite_gradient_theta_eta(1), diff_theta_eta(1, 0));

  Eigen::VectorXd W_root = (-hessian).cwiseSqrt();
  Eigen::MatrixXd diff2_theta_eta
    = diff_functor.diff2_theta_eta(theta, eta, W_root);

  Eigen::VectorXd finite_hessian_theta_eta
   = (hessian_theta_u - hessian_theta_l) / (2 * epsilon);

  EXPECT_FLOAT_EQ(finite_hessian_theta_eta(0), diff2_theta_eta(0, 0));
  EXPECT_FLOAT_EQ(finite_hessian_theta_eta(1), diff2_theta_eta(1, 0));
}

// unit tests for derivatives of B
template <typename T>
Eigen::MatrixXd compute_B(const Eigen::VectorXd& theta,
                          const Eigen::VectorXd& eta,
                          const Eigen::MatrixXd& covariance,
                          T diff_functor) {
  int group_size = theta.size();
  Eigen::VectorXd l_grad, hessian;
  diff_functor.diff(theta, eta, l_grad, hessian);
  Eigen::VectorXd W_root = (- hessian).cwiseSqrt();

  return Eigen::MatrixXd::Identity(group_size, group_size)
    + stan::math::quad_form_diag(covariance, W_root);
}

TEST(laplace, neg_binomial_2_log_dbl) {
  using stan::math::to_vector;
  using stan::math::diff_neg_binomial_2_log;
  using stan::math::sqr_exp_kernel_functor;
  using stan::math::laplace_marginal_density;
  using stan::math::var;
  using stan::math::value_of;

  int dim_phi = 2, dim_eta = 1, dim_theta = 2;
  Eigen::VectorXd phi(dim_phi), eta(dim_eta), theta_0(dim_theta);
  phi << 1.6, 0.45;
  eta << 1;
  theta_0 << 0, 0;

  std::vector<Eigen::VectorXd> x(dim_theta);
  Eigen::VectorXd x_0(2), x_1(2);
  x_0 <<  0.05100797, 0.16086164;
  x_1 << -0.59823393, 0.98701425;
  x[0] = x_0;
  x[1] = x_1;

  std::vector<double> delta;
  std::vector<int> delta_int;
  std::vector<int> y_index = {0, 1};
  Eigen::VectorXd y = to_vector({1, 6});

  diff_neg_binomial_2_log diff_functor(y, y_index, dim_theta);
  stan::math::sqr_exp_kernel_functor K;

  double log_p = laplace_marginal_density(diff_functor, K, phi, eta, x, delta,
                                          delta_int, theta_0);

  Eigen::Matrix<var, Eigen::Dynamic, 1> phi_v = phi, eta_v = eta;

  var target
    = laplace_marginal_density(diff_functor, K, phi_v, eta_v, x, delta,
                               delta_int, theta_0);

  VEC g;
  AVEC parm_vec = createAVEC(phi_v(0), phi_v(1), eta_v(0));
  target.grad(parm_vec, g);

  std::cout << "autodiff: ";
  for (size_t i = 0; i < g.size(); i++) std::cout << g[i] << " ";
  std::cout << std::endl;

  // finite diff test
  double diff = 1e-7;
  Eigen::VectorXd phi_dbl = value_of(phi), eta_dbl = value_of(eta);
  Eigen::VectorXd phi_1l = phi_dbl, phi_1u = phi_dbl,
    phi_2l = phi_dbl, phi_2u = phi_dbl, eta_l = eta_dbl, eta_u = eta_dbl;
  phi_1l(0) -= diff;
  phi_1u(0) += diff;
  phi_2l(1) -= diff;
  phi_2u(1) += diff;
  eta_l(0) -= diff;
  eta_u(0) += diff;

  double target_phi_1u = laplace_marginal_density(diff_functor, K, phi_1u,
                                                         eta_dbl, x, delta,
                                                         delta_int, theta_0),
         target_phi_1l = laplace_marginal_density(diff_functor, K, phi_1l,
                                                         eta_dbl, x, delta,
                                                         delta_int, theta_0),
         target_phi_2u = laplace_marginal_density(diff_functor, K, phi_2u,
                                                         eta_dbl, x, delta,
                                                         delta_int, theta_0),
         target_phi_2l = laplace_marginal_density(diff_functor, K, phi_2l,
                                                         eta_dbl, x, delta,
                                                         delta_int, theta_0),
         target_eta_u = laplace_marginal_density(diff_functor, K, phi_dbl,
                                                         eta_u, x, delta,
                                                         delta_int, theta_0),
         target_eta_l = laplace_marginal_density(diff_functor, K, phi_dbl,
                                                         eta_l, x, delta,
                                                         delta_int, theta_0);

  VEC g_finite(dim_phi + dim_eta);
  g_finite[0] = (target_phi_1u - target_phi_1l) / (2 * diff);
  g_finite[1] = (target_phi_2u - target_phi_2l) / (2 * diff);
  g_finite[2] = (target_eta_u - target_eta_l) / (2 * diff);

  std::cout << "Finite: ";
  for (int i = 0; i < 3; i++) std::cout << g_finite[i] << " ";
  std::cout << std::endl;

  // Save relevant variables for more detailed unit tests.
  Eigen::MatrixXd covariance, L;
  Eigen::VectorXd theta, a, theta_u, theta_l, W_root, l_grad, hessian;
  target =  laplace_marginal_density(diff_functor, K, phi, eta_dbl, x,
      delta, delta_int, covariance, theta, W_root, L, a, l_grad, theta_0);

  Eigen::MatrixXd B = compute_B(theta, eta_dbl, covariance, diff_functor),
  B_l = compute_B(theta, eta_l, covariance, diff_functor),
  B_u = compute_B(theta, eta_u, covariance, diff_functor);

  std::cout << std::endl << "B finite diff tests" << std::endl;
  std::cout << "log|B|: " << log(B.determinant()) << std::endl;
  std::cout << "finite diff: "
    << (log(B_u.determinant()) - log(B_l.determinant())) / (2 * diff)
            << std::endl;

  diff_functor.diff(theta, eta_dbl, l_grad, hessian);
  W_root = (-hessian).cwiseSqrt();

  Eigen::VectorXd hessian_l, hessian_u;
  diff_functor.diff(theta, eta_l, l_grad, hessian_l);
  diff_functor.diff(theta, eta_u, l_grad, hessian_u);

  // L = stan::math::cholesky_decompose(B);

  // std::cout << "candiate 1: " <<
  //   (B.inverse() * (-covariance * diff_functor.diff2_theta_eta(theta_0, eta_dbl, W_root)
  //    )).trace() << std::endl;

  Eigen::VectorXd W_finite_diff
    = (hessian_u - hessian_l) / (2 * diff);

  Eigen::VectorXd W_root_finite_diff
    = ((-hessian_u).cwiseSqrt() - (-hessian_l).cwiseSqrt()) / (2 * diff);

  Eigen::VectorXd W_root_diff = - 0.5 * stan::math::elt_divide(
    diff_functor.diff2_theta_eta(theta, eta_dbl, W_root), W_root);

  std::cout << "candiate 2: " <<
     2 * (B.inverse() * (W_root.asDiagonal() * covariance * W_root_diff.asDiagonal())
   ).trace() << std::endl;

  double diff_log_B =
    - (L.transpose().triangularView<Eigen::Upper>()
      .solve(L.triangularView<Eigen::Lower>()
       .solve(W_root.asDiagonal() * covariance
      * stan::math::elt_divide(diff_functor.
        diff2_theta_eta(theta, eta_dbl, W_root), W_root).asDiagonal()))).trace();

  std::cout << "candiate 3: " << diff_log_B << std::endl;
  std::cout << "full log B term: " << - 0.5 * diff_log_B << std::endl;

  std::cout << std::endl << "W_root: " << W_root.transpose() << std::endl;

  std::cout << "W_root finite diff: " << W_root_finite_diff.transpose() << std::endl;
  std::cout << "W_root diff: " << W_root_diff.transpose() << std::endl;

  std::cout << std::endl << "diff2 finite: " << W_finite_diff.transpose() << std::endl;
  std::cout << std::endl << "diff: " <<
    diff_functor.diff2_theta_eta(theta, eta_dbl, W_root).transpose() << std::endl;

  std::cout << std::endl << "Differentiation of theta star." << std::endl;


  Eigen::VectorXd b = covariance * diff_functor.diff_theta_eta(theta, eta_dbl);
    Eigen::VectorXd s3 = (Eigen::MatrixXd::Identity(theta.size(), theta.size())
      + covariance * stan::math::square(W_root).asDiagonal()).inverse() * b;

  Eigen::VectorXd theta_star = theta;

  target = laplace_marginal_density(diff_functor, K, phi, eta_u, x,
      delta, delta_int, covariance, theta_u, W_root, L, a, l_grad, theta_0);
  target = laplace_marginal_density(diff_functor, K, phi, eta_l, x,
      delta, delta_int, covariance, theta_l, W_root, L, a, l_grad, theta_0);

  std::cout << "theta: " << theta.transpose() << std::endl;
  std::cout << "theta finite diff: "
    << ((theta_u - theta_l) / (2 * diff)).transpose()
    << std::endl;

  std::cout << "theta analytical diff: " << s3.transpose() << std::endl;

  std::cout << std::endl << "Computation of s2: " << std::endl;

  // Reset the variables to their unperturbed states.
  target =  laplace_marginal_density(diff_functor, K, phi, eta_dbl, x,
      delta, delta_int, covariance, theta, W_root, L, a, l_grad, theta_0);

  Eigen::VectorXd theta_1l = theta_star, theta_1u = theta_star,
                  theta_2l = theta_star, theta_2u = theta_star;
  theta_1l(0) -= diff;
  theta_1u(0) += diff;
  theta_2l(1) -= diff;
  theta_2u(1) += diff;

  B = compute_B(theta_star, eta, covariance, diff_functor);

  double
    log_B_1l = log(compute_B(theta_1l, eta_dbl, covariance, diff_functor).determinant()),
    log_B_1u = log(compute_B(theta_1u, eta_dbl, covariance, diff_functor).determinant()),
    log_B_2l = log(compute_B(theta_2l, eta_dbl, covariance, diff_functor).determinant()),
    log_B_2u = log(compute_B(theta_2u, eta_dbl, covariance, diff_functor).determinant());

  Eigen::VectorXd s2_finite(2);
  s2_finite(0) = - 0.5 * (log_B_1u - log_B_1l) / (2 * diff);
  s2_finite(1) = - 0.5 * (log_B_2u - log_B_2l) / (2 * diff);

  std::cout << "s2_finite: " << s2_finite.transpose() << std::endl;


  std::cout << "log diff finite: " <<
    (diff_functor.log_likelihood(theta_star, eta_u)
      - diff_functor.log_likelihood(theta_star, eta_l)) / (2 * diff)
      << std::endl;


  // double tol = 1e-4;
  // EXPECT_NEAR(g_finite[0], g[0], tol);
  // EXPECT_NEAR(g_finite[1], g[1], tol);
  // EXPECT_NEAR(g_finite[2], g[2], tol);
}
