#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <test/unit/util.hpp>
#include <test/unit/math/rev/prob/expect_eq_diffs.hpp>
#include <test/unit/math/rev/prob/test_gradients.hpp>
#include <test/unit/math/rev/prob/test_gradients_multi_student_t_cholesky.hpp>
#include <test/unit/math/prim/prob/agrad_distributions_multi_student_t_cholesky.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <string>

template <typename T_y, typename T_dof, typename T_loc, typename T_scale>
void expect_propto_multi_student_t_cholesky_lpdf(T_y y1, T_dof nu1, T_loc mu1,
                                                 T_scale L1, T_y y2, T_dof nu2,
                                                 T_loc mu2, T_scale L2,
                                                 std::string message = "") {
  expect_eq_diffs(
      stan::math::multi_student_t_cholesky_lpdf<false>(y1, nu1, mu1, L1),
      stan::math::multi_student_t_cholesky_lpdf<false>(y2, nu2, mu2, L2),
      stan::math::multi_student_t_cholesky_lpdf<true>(y1, nu1, mu1, L1),
      stan::math::multi_student_t_cholesky_lpdf<true>(y2, nu2, mu2, L2),
      message);
}

TEST_F(agrad_distributions_multi_student_t_cholesky, Propto) {
  using stan::math::to_var;
  expect_propto_multi_student_t_cholesky_lpdf(
      to_var(y), to_var(nu), to_var(mu), to_var(L), to_var(y2), to_var(nu),
      to_var(mu2), to_var(L2), "All vars: y, nu, mu, L");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_student_t_cholesky, ProptoY) {
  using stan::math::to_var;
  expect_propto_multi_student_t_cholesky_lpdf(to_var(y), nu, mu, L, to_var(y2),
                                              nu, mu, L, "var: y");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_student_t_cholesky, ProptoYMu) {
  using stan::math::to_var;
  expect_propto_multi_student_t_cholesky_lpdf(to_var(y), nu, to_var(mu), L,
                                              to_var(y2), nu, to_var(mu2), L,
                                              "var: y and mu");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_student_t_cholesky, ProptoYL) {
  using stan::math::to_var;
  expect_propto_multi_student_t_cholesky_lpdf(to_var(y), nu, mu, to_var(L),
                                              to_var(y2), nu, mu, to_var(L2),
                                              "var: y and L");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_student_t_cholesky, ProptoMu) {
  using stan::math::to_var;
  expect_propto_multi_student_t_cholesky_lpdf(y, nu, to_var(mu), L, y, nu,
                                              to_var(mu2), L, "var: mu");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_student_t_cholesky, ProptoMuL) {
  using stan::math::to_var;
  expect_propto_multi_student_t_cholesky_lpdf(y, nu, to_var(mu), to_var(L), y,
                                              nu, to_var(mu2), to_var(L2),
                                              "var: mu and L");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_student_t_cholesky, ProptoL) {
  using stan::math::to_var;
  expect_propto_multi_student_t_cholesky_lpdf(y, nu, mu, to_var(L), y, nu, mu,
                                              to_var(L2), "var: L");

  stan::math::recover_memory();
}

TEST(ProbDistributionsMultiStudentTCholesky, MultiStudentTVar) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::math::var;
  using std::vector;
  var nu(5);
  Matrix<var, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<var, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<var, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  Matrix<var, Dynamic, Dynamic> L = Sigma.llt().matrixL();

  EXPECT_FLOAT_EQ(
      -10.2136949646,
      stan::math::multi_student_t_cholesky_lpdf(y, nu, mu, L).val());

  stan::math::recover_memory();
}
TEST(ProbDistributionsMultiStudentTCholesky, MultiStudentTGradientUnivariate) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using Eigen::VectorXd;
  using stan::math::multi_student_t_cholesky_lpdf;
  using stan::math::var;
  using std::vector;

  var nu_var(5);

  Matrix<var, Dynamic, 1> y_var(1, 1);
  y_var << 2.0;

  Matrix<var, Dynamic, 1> mu_var(1, 1);
  mu_var << 1.0;

  Matrix<var, Dynamic, Dynamic> L_var(1, 1);
  L_var(0, 0) = 3.0;

  std::vector<var> x;
  x.push_back(y_var(0));
  x.push_back(mu_var(0));
  x.push_back(L_var(0, 0));
  x.push_back(nu_var);

  var lp
      = stan::math::multi_student_t_cholesky_lpdf(y_var, nu_var, mu_var, L_var);
  vector<double> grad;
  lp.grad(x, grad);

  // ===================================

  double nu(5);

  Matrix<double, Dynamic, 1> y(1, 1);
  y << 2.0;

  Matrix<double, Dynamic, 1> mu(1, 1);
  mu << 1.0;

  Matrix<double, Dynamic, Dynamic> L(1, 1);
  L << 3.0;

  double epsilon = 1e-6;

  Matrix<double, Dynamic, 1> y_m(1, 1);
  Matrix<double, Dynamic, 1> y_p(1, 1);
  y_p[0] = y[0] + epsilon;
  y_m[0] = y[0] - epsilon;
  double grad_diff = (multi_student_t_cholesky_lpdf(y_p, nu, mu, L)
                      - multi_student_t_cholesky_lpdf(y_m, nu, mu, L))
                     / (2 * epsilon);
  EXPECT_FLOAT_EQ(grad_diff, grad[0]);

  Matrix<double, Dynamic, 1> mu_m(1, 1);
  Matrix<double, Dynamic, 1> mu_p(1, 1);
  mu_p[0] = mu[0] + epsilon;
  mu_m[0] = mu[0] - epsilon;
  grad_diff = (multi_student_t_cholesky_lpdf(y, nu, mu_p, L)
               - multi_student_t_cholesky_lpdf(y, nu, mu_m, L))
              / (2 * epsilon);
  EXPECT_FLOAT_EQ(grad_diff, grad[1]);

  Matrix<double, Dynamic, Dynamic> L_m(1, 1);
  Matrix<double, Dynamic, Dynamic> L_p(1, 1);
  L_p(0) = L(0) + epsilon;
  L_m(0) = L(0) - epsilon;
  grad_diff = (multi_student_t_cholesky_lpdf(y, nu, mu, L_p)
               - multi_student_t_cholesky_lpdf(y, nu, mu, L_m))
              / (2 * epsilon);
  EXPECT_FLOAT_EQ(grad_diff, grad[2]);

  double nu_p(nu + epsilon);
  double nu_m(nu - epsilon);
  grad_diff = (multi_student_t_cholesky_lpdf(y, nu_p, mu, L)
               - multi_student_t_cholesky_lpdf(y, nu_m, mu, L))
              / (2 * epsilon);
  EXPECT_FLOAT_EQ(grad_diff, grad[3]);

  stan::math::recover_memory();
}

struct multi_student_t_cholesky_fun {
  const int K_;

  explicit multi_student_t_cholesky_fun(int K) : K_(K) {}

  template <typename T>
  T operator()(const std::vector<T>& x) const {
    using Eigen::Dynamic;
    using Eigen::Matrix;
    using stan::math::var;
    T nu;
    Matrix<T, Dynamic, 1> y(K_);
    Matrix<T, Dynamic, 1> mu(K_);
    Matrix<T, Dynamic, Dynamic> Sigma(K_, K_);
    int pos = 0;
    for (int i = 0; i < K_; ++i)
      y(i) = x[pos++];
    for (int i = 0; i < K_; ++i)
      mu(i) = x[pos++];
    for (int j = 0; j < K_; ++j) {
      for (int i = 0; i <= j; ++i) {
        Sigma(i, j) = x[pos++];
        Sigma(j, i) = Sigma(i, j);
      }
    }
    Matrix<T, Dynamic, Dynamic> L = Sigma.llt().matrixL();
    nu = x[pos++];
    return stan::math::multi_student_t_cholesky_lpdf<false>(y, nu, mu, L);
  }
};

TEST(ProbDistributionsMultiStudentTCholesky, TestGradFunctional) {
  std::vector<double> x(3 + 3 + 3 * 2 + 1);
  // y
  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = -3.0;
  // mu
  x[3] = 0.0;
  x[4] = -2.0;
  x[5] = -3.0;
  // Sigma
  x[6] = 1;
  x[7] = -1;
  x[8] = 10;
  x[9] = -2;
  x[10] = 20;
  x[11] = 56;
  // nu
  x[12] = 5;

  test_grad(multi_student_t_cholesky_fun(3), x);

  std::vector<double> u(4);
  u[0] = 1.9;
  u[1] = -2.7;
  u[2] = 0.48;
  u[3] = 5;

  test_grad(multi_student_t_cholesky_fun(1), u);

  stan::math::recover_memory();
}

template <int is_row_vec_y, int is_row_vec_mu>
struct vectorized_multi_student_t_cholesky_fun {
  // size of each vector and order of square matrix sigma
  // size of the array of eigen vectors
  // direct use eigen vector for y
  // direct use eigen vector for mu
  const int K_;
  const int L_;
  const bool dont_vectorize_y;
  const bool dont_vectorize_mu;

  vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(
      int K, int L, bool M = false, bool N = false)
      : K_(K), L_(L), dont_vectorize_y(M), dont_vectorize_mu(N) {
    if ((dont_vectorize_y || dont_vectorize_mu) && L != 1)
      throw std::runtime_error(
          "attempt to disable vectorization with vector "
          "bigger than 1");
  }

  template <typename T_y, typename T_mu, typename T_sigma, typename T_nu>
  stan::return_type_t<T_y, T_mu, T_sigma, T_nu> operator()(
      const std::vector<T_y>& y_vec, const std::vector<T_mu>& mu_vec,
      const std::vector<T_sigma>& sigma_vec, const T_nu& nu) const {
    using Eigen::Dynamic;
    using Eigen::Matrix;
    using std::vector;
    vector<Matrix<T_y, is_row_vec_y, is_row_vec_y * -1> > y(
        L_, Matrix<T_y, is_row_vec_y, is_row_vec_y * -1>(K_));
    vector<Matrix<T_mu, is_row_vec_mu, is_row_vec_mu * -1> > mu(
        L_, Matrix<T_mu, is_row_vec_mu, is_row_vec_mu * -1>(K_));
    Matrix<T_sigma, Dynamic, Dynamic> Sigma(K_, K_);
    int pos = 0;
    for (int i = 0; i < L_; ++i)
      for (int j = 0; j < K_; ++j)
        y[i](j) = y_vec[pos++];

    pos = 0;
    for (int i = 0; i < L_; ++i)
      for (int j = 0; j < K_; ++j)
        mu[i](j) = mu_vec[pos++];

    pos = 0;
    for (int j = 0; j < K_; ++j) {
      for (int i = 0; i <= j; ++i) {
        Sigma(i, j) = sigma_vec[pos++];
        Sigma(j, i) = Sigma(i, j);
      }
    }
    Matrix<T_sigma, Dynamic, Dynamic> L = Sigma.llt().matrixL();

    if (dont_vectorize_y) {
      if (dont_vectorize_mu)
        return stan::math::multi_student_t_cholesky_lpdf<false>(y[0], nu, mu[0],
                                                                L);
      else
        return stan::math::multi_student_t_cholesky_lpdf<false>(y[0], nu, mu,
                                                                L);
    } else {
      if (dont_vectorize_mu)
        return stan::math::multi_student_t_cholesky_lpdf<false>(y, nu, mu[0],
                                                                L);
      else
        return stan::math::multi_student_t_cholesky_lpdf<false>(y, nu, mu, L);
    }
  }
};

template <int is_row_vec_y, int is_row_vec_mu>
void test_all_multi_student_t_cholesky() {
  {
    using stan::math::var;
    std::vector<double> y_(3), mu_(3), sigma_(6);
    // y
    y_[0] = 1.0;
    y_[1] = 2.0;
    y_[2] = -3.0;
    // mu
    mu_[0] = 0.0;
    mu_[1] = -2.0;
    mu_[2] = -3.0;
    // Sigma
    sigma_[0] = 1;
    sigma_[1] = -1;
    sigma_[2] = 10;
    sigma_[3] = -2;
    sigma_[4] = 20;
    sigma_[5] = 56;
    for (int ii = 0; ii < 2; ii++)
      for (int jj = 0; jj < 2; jj++) {
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            y_, mu_, sigma_, 5);
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            y_, mu_, stan::math::to_var(sigma_), 5);
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            y_, stan::math::to_var(mu_), sigma_, 5);
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            y_, stan::math::to_var(mu_), stan::math::to_var(sigma_), 5);
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            stan::math::to_var(y_), mu_, sigma_, 5);
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            stan::math::to_var(y_), mu_, stan::math::to_var(sigma_), 5);
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            stan::math::to_var(y_), stan::math::to_var(mu_), sigma_, 5);
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            stan::math::to_var(y_), stan::math::to_var(mu_),
            stan::math::to_var(sigma_), 5);
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            y_, mu_, sigma_, var(5));
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            y_, mu_, stan::math::to_var(sigma_), var(5));
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            y_, stan::math::to_var(mu_), sigma_, var(5));
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            y_, stan::math::to_var(mu_), stan::math::to_var(sigma_), var(5));
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            stan::math::to_var(y_), mu_, sigma_, var(5));
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            stan::math::to_var(y_), mu_, stan::math::to_var(sigma_), var(5));
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            stan::math::to_var(y_), stan::math::to_var(mu_), sigma_, var(5));
        test_grad_multi_student_t_cholesky(
            vectorized_multi_student_t_cholesky_fun<is_row_vec_y,
                                                    is_row_vec_mu>(3, 1, ii,
                                                                   jj),
            stan::math::to_var(y_), stan::math::to_var(mu_),
            stan::math::to_var(sigma_), var(5));
      }
  }

  {
    using stan::math::var;
    std::vector<double> y_(6), mu_(6), sigma_(6);
    // y[1]
    y_[0] = 1.0;
    y_[1] = 2.0;
    y_[2] = -3.0;
    // y[2]
    y_[3] = 0.0;
    y_[4] = -2.0;
    y_[5] = -3.0;

    // mu[1]
    mu_[0] = 0.0;
    mu_[1] = 1.0;
    mu_[2] = 3.0;
    // mu[2]
    mu_[3] = 0.0;
    mu_[4] = -1.0;
    mu_[5] = -2.0;

    // Sigma
    sigma_[0] = 1;
    sigma_[1] = -1;
    sigma_[2] = 10;
    sigma_[3] = -2;
    sigma_[4] = 20;
    sigma_[5] = 56;

    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        y_, mu_, sigma_, 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        y_, mu_, stan::math::to_var(sigma_), 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        y_, stan::math::to_var(mu_), sigma_, 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        y_, stan::math::to_var(mu_), stan::math::to_var(sigma_), 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        stan::math::to_var(y_), mu_, sigma_, 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        stan::math::to_var(y_), mu_, stan::math::to_var(sigma_), 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        stan::math::to_var(y_), stan::math::to_var(mu_), sigma_, 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        stan::math::to_var(y_), stan::math::to_var(mu_),
        stan::math::to_var(sigma_), 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        y_, mu_, sigma_, var(5));
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        y_, mu_, stan::math::to_var(sigma_), var(5));
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        y_, stan::math::to_var(mu_), sigma_, var(5));
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        y_, stan::math::to_var(mu_), stan::math::to_var(sigma_), var(5));
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        stan::math::to_var(y_), mu_, sigma_, var(5));
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        stan::math::to_var(y_), mu_, stan::math::to_var(sigma_), var(5));
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        stan::math::to_var(y_), stan::math::to_var(mu_), sigma_, var(5));
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(3,
                                                                             2),
        stan::math::to_var(y_), stan::math::to_var(mu_),
        stan::math::to_var(sigma_), var(5));
  }
  {
    using stan::math::var;
    std::vector<double> y_(1), mu_(1), sigma_(1);
    y_[0] = 1.9;
    mu_[0] = -2.7;
    sigma_[0] = 0.48;

    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        y_, mu_, sigma_, 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        y_, mu_, stan::math::to_var(sigma_), 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        y_, stan::math::to_var(mu_), sigma_, 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        y_, stan::math::to_var(mu_), stan::math::to_var(sigma_), 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        stan::math::to_var(y_), mu_, sigma_, 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        stan::math::to_var(y_), mu_, stan::math::to_var(sigma_), 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        stan::math::to_var(y_), stan::math::to_var(mu_), sigma_, 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        stan::math::to_var(y_), stan::math::to_var(mu_),
        stan::math::to_var(sigma_), 5);
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        y_, mu_, sigma_, var(5));
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        y_, mu_, stan::math::to_var(sigma_), var(5));
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        y_, stan::math::to_var(mu_), sigma_, var(5));
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        y_, stan::math::to_var(mu_), stan::math::to_var(sigma_), var(5));
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        stan::math::to_var(y_), mu_, sigma_, var(5));
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        stan::math::to_var(y_), mu_, stan::math::to_var(sigma_), var(5));
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        stan::math::to_var(y_), stan::math::to_var(mu_), sigma_, var(5));
    test_grad_multi_student_t_cholesky(
        vectorized_multi_student_t_cholesky_fun<is_row_vec_y, is_row_vec_mu>(1,
                                                                             1),
        stan::math::to_var(y_), stan::math::to_var(mu_),
        stan::math::to_var(sigma_), var(5));
  }
}

TEST(ProbDistributionsMultiStudentTCholesky, TestGradFunctionalVectorized) {
  test_all_multi_student_t_cholesky<1, 1>();
  test_all_multi_student_t_cholesky<1, -1>();
  test_all_multi_student_t_cholesky<-1, 1>();
  test_all_multi_student_t_cholesky<-1, -1>();

  stan::math::recover_memory();
}

TEST(ProbDistributionsMultiStudentTCholesky, check_varis_on_stack) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using std::vector;
  double nu(5);
  Matrix<double, Dynamic, 1> y(3, 1);
  y << 2.0, -2.0, 11.0;
  Matrix<double, Dynamic, 1> mu(3, 1);
  mu << 1.0, -1.0, 3.0;
  Matrix<double, Dynamic, Dynamic> Sigma(3, 3);
  Sigma << 9.0, -3.0, 0.0, -3.0, 4.0, 0.0, 0.0, 0.0, 5.0;
  Matrix<double, Dynamic, Dynamic> L = Sigma.llt().matrixL();

  using stan::math::multi_student_t_cholesky_lpdf;
  using stan::math::to_var;
  test::check_varis_on_stack(multi_student_t_cholesky_lpdf<false>(
      to_var(y), to_var(nu), to_var(mu), to_var(L)));
  test::check_varis_on_stack(multi_student_t_cholesky_lpdf<false>(
      to_var(y), to_var(nu), to_var(mu), L));
  test::check_varis_on_stack(multi_student_t_cholesky_lpdf<false>(
      to_var(y), to_var(nu), mu, to_var(L)));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<false>(to_var(y), to_var(nu), mu, L));
  test::check_varis_on_stack(multi_student_t_cholesky_lpdf<false>(
      to_var(y), nu, to_var(mu), to_var(L)));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<false>(to_var(y), nu, to_var(mu), L));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<false>(to_var(y), nu, mu, to_var(L)));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<false>(to_var(y), nu, mu, L));
  test::check_varis_on_stack(multi_student_t_cholesky_lpdf<false>(
      y, to_var(nu), to_var(mu), to_var(L)));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<false>(y, to_var(nu), to_var(mu), L));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<false>(y, to_var(nu), mu, to_var(L)));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<false>(y, to_var(nu), mu, L));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<false>(y, nu, to_var(mu), to_var(L)));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<false>(y, nu, to_var(mu), L));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<false>(y, nu, mu, to_var(L)));
  test::check_varis_on_stack(multi_student_t_cholesky_lpdf<true>(
      to_var(y), to_var(nu), to_var(mu), to_var(L)));
  test::check_varis_on_stack(multi_student_t_cholesky_lpdf<true>(
      to_var(y), to_var(nu), to_var(mu), L));
  test::check_varis_on_stack(multi_student_t_cholesky_lpdf<true>(
      to_var(y), to_var(nu), mu, to_var(L)));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<true>(to_var(y), to_var(nu), mu, L));
  test::check_varis_on_stack(multi_student_t_cholesky_lpdf<true>(
      to_var(y), nu, to_var(mu), to_var(L)));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<true>(to_var(y), nu, to_var(mu), L));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<true>(to_var(y), nu, mu, to_var(L)));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<true>(to_var(y), nu, mu, L));
  test::check_varis_on_stack(multi_student_t_cholesky_lpdf<true>(
      y, to_var(nu), to_var(mu), to_var(L)));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<true>(y, to_var(nu), to_var(mu), L));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<true>(y, to_var(nu), mu, to_var(L)));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<true>(y, to_var(nu), mu, L));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<true>(y, nu, to_var(mu), to_var(L)));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<true>(y, nu, to_var(mu), L));
  test::check_varis_on_stack(
      multi_student_t_cholesky_lpdf<true>(y, nu, mu, to_var(L)));

  stan::math::recover_memory();
}
