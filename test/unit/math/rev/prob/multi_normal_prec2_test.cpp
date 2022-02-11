#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/prob/test_gradients.hpp>
#include <test/unit/math/rev/prob/test_gradients_multi_normal.hpp>
#include <test/unit/math/rev/prob/expect_eq_diffs.hpp>
#include <test/unit/math/prim/prob/agrad_distributions_multi_normal_multi_row.hpp>
#include <test/unit/math/prim/prob/agrad_distributions_multi_normal.hpp>
#include <vector>
#include <string>

template <typename T_y, typename T_loc, typename T_scale>
void expect_propto_multi_normal_prec_log(T_y y1, T_loc mu1, T_scale sigma1,
                                         T_y y2, T_loc mu2, T_scale sigma2,
                                         std::string message = "") {
  expect_eq_diffs(stan::math::multi_normal_prec_log<false>(y1, mu1, sigma1),
                  stan::math::multi_normal_prec_log<false>(y2, mu2, sigma2),
                  stan::math::multi_normal_prec_log<true>(y1, mu1, sigma1),
                  stan::math::multi_normal_prec_log<true>(y2, mu2, sigma2),
                  message);
}

TEST_F(agrad_distributions_multi_normal, Propto_prec2) {
  using stan::math::to_var;
  expect_propto_multi_normal_prec_log(to_var(y), to_var(mu), to_var(Sigma),
                                      to_var(y2), to_var(mu2), to_var(Sigma2),
                                      "All vars: y, mu, sigma");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_normal, ProptoY_prec2) {
  using stan::math::to_var;
  expect_propto_multi_normal_prec_log(to_var(y), mu, Sigma, to_var(y2), mu,
                                      Sigma, "var: y");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_normal, ProptoYMu_prec2) {
  using stan::math::to_var;
  expect_propto_multi_normal_prec_log(to_var(y), to_var(mu), Sigma, to_var(y2),
                                      to_var(mu2), Sigma, "var: y and mu");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_normal, ProptoYSigma_prec2) {
  using stan::math::to_var;
  expect_propto_multi_normal_prec_log(to_var(y), mu, to_var(Sigma), to_var(y2),
                                      mu, to_var(Sigma2), "var: y and sigma");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_normal, ProptoMu_prec2) {
  using stan::math::to_var;
  expect_propto_multi_normal_prec_log(y, to_var(mu), Sigma, y, to_var(mu2),
                                      Sigma, "var: mu");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_normal, ProptoMuSigma_prec2) {
  using stan::math::to_var;
  expect_propto_multi_normal_prec_log(y, to_var(mu), to_var(Sigma), y,
                                      to_var(mu2), to_var(Sigma2),
                                      "var: mu and sigma");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_normal, ProptoSigma_prec2) {
  using stan::math::to_var;
  expect_propto_multi_normal_prec_log(y, mu, to_var(Sigma), y, mu,
                                      to_var(Sigma2), "var: sigma");

  stan::math::recover_memory();
}

TEST_F(agrad_distributions_multi_normal_multi_row, Propto_prec2) {
  using stan::math::to_var;
  expect_propto_multi_normal_prec_log(to_var(y), to_var(mu), to_var(Sigma),
                                      to_var(y2), to_var(mu2), to_var(Sigma2),
                                      "All vars: y, mu, sigma");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_normal_multi_row, ProptoY_prec2) {
  using stan::math::to_var;
  expect_propto_multi_normal_prec_log(to_var(y), mu, Sigma, to_var(y2), mu,
                                      Sigma, "var: y");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_normal_multi_row, ProptoYMu_prec2) {
  using stan::math::to_var;
  expect_propto_multi_normal_prec_log(to_var(y), to_var(mu), Sigma, to_var(y2),
                                      to_var(mu2), Sigma, "var: y and mu");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_normal_multi_row, ProptoYSigma_prec2) {
  using stan::math::to_var;
  expect_propto_multi_normal_prec_log(to_var(y), mu, to_var(Sigma), to_var(y2),
                                      mu, to_var(Sigma2), "var: y and sigma");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_normal_multi_row, ProptoMu_prec2) {
  using stan::math::to_var;
  expect_propto_multi_normal_prec_log(y, to_var(mu), Sigma, y, to_var(mu2),
                                      Sigma, "var: mu");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_normal_multi_row, ProptoMuSigma_prec2) {
  using stan::math::to_var;
  expect_propto_multi_normal_prec_log(y, to_var(mu), to_var(Sigma), y,
                                      to_var(mu2), to_var(Sigma2),
                                      "var: mu and sigma");

  stan::math::recover_memory();
}
TEST_F(agrad_distributions_multi_normal_multi_row, ProptoSigma_prec2) {
  using stan::math::to_var;
  expect_propto_multi_normal_prec_log(y, mu, to_var(Sigma), y, mu,
                                      to_var(Sigma2), "var: sigma");

  stan::math::recover_memory();
}

struct multi_normal_prec_fun {
  const int K_;

  explicit multi_normal_prec_fun(int K) : K_(K) {}

  template <typename T>
  T operator()(const std::vector<T>& x) const {
    using Eigen::Dynamic;
    using Eigen::Matrix;
    using stan::math::var;
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
    return stan::math::multi_normal_prec_log<false>(y, mu, Sigma);
  }
};

TEST(ProbDistributionsMultiNormalPrec2, TestGradFunctional) {
  std::vector<double> x(3 + 3 + 3 * 2);
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

  test_grad(multi_normal_prec_fun(3), x);

  std::vector<double> u(3);
  u[0] = 1.9;
  u[1] = -2.7;
  u[2] = 0.48;

  test_grad(multi_normal_prec_fun(1), u);

  stan::math::recover_memory();
}

template <int is_row_vec_y, int is_row_vec_mu>
struct vectorized_multi_normal_prec_fun {
  // size of each vector and order of square matrix sigma
  const int K_;
  // size of the array of eigen vectors
  const int L_;
  // direct use eigen vector for y
  const bool dont_vectorize_y;
  // direct use eigen vector for mu
  const bool dont_vectorize_mu;

  vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(int K, int L,
                                                                bool M = false,
                                                                bool N = false)
      : K_(K), L_(L), dont_vectorize_y(M), dont_vectorize_mu(N) {
    if ((dont_vectorize_y || dont_vectorize_mu) && L != 1)
      throw std::runtime_error(
          "attempt to disable vectorization with vector "
          "bigger than 1");
  }

  template <typename T_y, typename T_mu, typename T_sigma>
  stan::return_type_t<T_y, T_mu, T_sigma> operator()(
      const std::vector<T_y>& y_vec, const std::vector<T_mu>& mu_vec,
      const std::vector<T_sigma>& sigma_vec) const {
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

    if (dont_vectorize_y) {
      if (dont_vectorize_mu)
        return stan::math::multi_normal_prec_log<false>(y[0], mu[0], Sigma);
      else
        return stan::math::multi_normal_prec_log<false>(y[0], mu, Sigma);
    } else {
      if (dont_vectorize_mu)
        return stan::math::multi_normal_prec_log<false>(y, mu[0], Sigma);
      else
        return stan::math::multi_normal_prec_log<false>(y, mu, Sigma);
    }
  }
};

template <int is_row_vec_y, int is_row_vec_mu>
void test_all_multi_normal_prec2() {
  {
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
        test_grad_multi_normal(
            vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(
                3, 1, ii, jj),
            y_, mu_, sigma_);
        test_grad_multi_normal(
            vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(
                3, 1, ii, jj),
            y_, mu_, stan::math::to_var(sigma_));
        test_grad_multi_normal(
            vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(
                3, 1, ii, jj),
            y_, stan::math::to_var(mu_), sigma_);
        test_grad_multi_normal(
            vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(
                3, 1, ii, jj),
            y_, stan::math::to_var(mu_), stan::math::to_var(sigma_));
        test_grad_multi_normal(
            vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(
                3, 1, ii, jj),
            stan::math::to_var(y_), mu_, sigma_);
        test_grad_multi_normal(
            vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(
                3, 1, ii, jj),
            stan::math::to_var(y_), mu_, stan::math::to_var(sigma_));
        test_grad_multi_normal(
            vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(
                3, 1, ii, jj),
            stan::math::to_var(y_), stan::math::to_var(mu_), sigma_);
        test_grad_multi_normal(
            vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(
                3, 1, ii, jj),
            stan::math::to_var(y_), stan::math::to_var(mu_),
            stan::math::to_var(sigma_));
      }
  }

  {
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

    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(3, 2), y_,
        mu_, sigma_);
    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(3, 2), y_,
        mu_, stan::math::to_var(sigma_));
    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(3, 2), y_,
        stan::math::to_var(mu_), sigma_);
    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(3, 2), y_,
        stan::math::to_var(mu_), stan::math::to_var(sigma_));
    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(3, 2),
        stan::math::to_var(y_), mu_, sigma_);
    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(3, 2),
        stan::math::to_var(y_), mu_, stan::math::to_var(sigma_));
    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(3, 2),
        stan::math::to_var(y_), stan::math::to_var(mu_), sigma_);
    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(3, 2),
        stan::math::to_var(y_), stan::math::to_var(mu_),
        stan::math::to_var(sigma_));
  }
  {
    std::vector<double> y_(1), mu_(1), sigma_(1);
    y_[0] = 1.9;
    mu_[0] = -2.7;
    sigma_[0] = 0.48;

    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(1, 1), y_,
        mu_, sigma_);
    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(1, 1), y_,
        mu_, stan::math::to_var(sigma_));
    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(1, 1), y_,
        stan::math::to_var(mu_), sigma_);
    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(1, 1), y_,
        stan::math::to_var(mu_), stan::math::to_var(sigma_));
    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(1, 1),
        stan::math::to_var(y_), mu_, sigma_);
    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(1, 1),
        stan::math::to_var(y_), mu_, stan::math::to_var(sigma_));
    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(1, 1),
        stan::math::to_var(y_), stan::math::to_var(mu_), sigma_);
    test_grad_multi_normal(
        vectorized_multi_normal_prec_fun<is_row_vec_y, is_row_vec_mu>(1, 1),
        stan::math::to_var(y_), stan::math::to_var(mu_),
        stan::math::to_var(sigma_));
  }
}

TEST(ProbDistributionsMultiNormalPrec2, TestGradFunctionalVectorized) {
  test_all_multi_normal_prec2<1, 1>();
  test_all_multi_normal_prec2<1, -1>();
  test_all_multi_normal_prec2<-1, 1>();
  test_all_multi_normal_prec2<-1, -1>();

  stan::math::recover_memory();
}
