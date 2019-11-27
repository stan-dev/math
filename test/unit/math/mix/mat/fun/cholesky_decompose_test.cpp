#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <vector>

// can't autodiff directly through Cholesky due to symmetry test;
// use unconstrained input and constrain to test Cholesky derivs;
// dof must be (n choose 2) + n
auto f(int dof) {
  return [=](const auto& x) {
    auto y = stan::math::cov_matrix_constrain(x, dof);
    return stan::math::cholesky_decompose(y);
  };
}

void expect_cholesky(const Eigen::MatrixXd& Sigma) {
  Eigen::VectorXd yy = stan::math::cov_matrix_free(Sigma);
  stan::test::expect_ad(f(10), yy);
}

TEST(MathMixMatFun, choleskyDecompose) {
  // 1 x 1 matrix;  (1 choose 2) + 1 = 1
  Eigen::VectorXd x1(1);
  x1 << 1;
  stan::test::expect_ad(f(1), x1);

  // 2 x 2 matrix;  (2 choose 2) + 2 = 3
  Eigen::VectorXd x3(3);
  x3 << 1, 2, -1;
  stan::test::expect_ad(f(2), x3);

  // 3 x 3 matrix;  (3 choose 2) + 3 = 6
  Eigen::VectorXd x6(6);
  x6 << 1, -1, 1.1, 1.4, 2.1, 0.7;
  stan::test::expect_ad(f(3), x6);

  // 4 x 4 matrix;  (4 choose 2) + 4 = 10
  Eigen::VectorXd x10(10);
  x10 << 1, -0.1, 1.1, 1.4, -1.1, 0.7, 1.0, 1.3, -0.5, 0.3;
  stan::test::expect_ad(f(4), x10);

  // 2 x 3 matrix will throw; test directly
  auto g = [](const auto& x) { return stan::math::cholesky_decompose(x); };
  Eigen::MatrixXd y(2, 3);
  y << 1, 2, 3, 4, 5, 6;
  stan::test::expect_ad(g, y);

  // asymmetric will throw
  Eigen::MatrixXd z(2, 2);
  z << 1, 2, 3, 4;
  stan::test::expect_ad(g, z);

  // general sizes
  for (int n = 0; n < 5; ++n) {
    int dof = n + (n * (n - 1)) / 2;
    Eigen::VectorXd y(dof);
    for (int i = 0; i < dof; ++i)
      y(i) = (i * 10) / 100.0;
    stan::test::expect_ad(f(dof), y);
  }

  // GP covar
  for (size_t n = 1; n < 5; ++n) {
    std::vector<double> xx(n);
    for (size_t i = 0; i < n; ++i) {
      xx[i] = (i * 10) / 100.0;
    }
    double alpha = 0.75;
    double length_scale = 1.25;
    double jitter = 0.1;
    Eigen::MatrixXd Sigma = stan::math::add_diag(
        stan::math::cov_exp_quad(xx, alpha, length_scale), jitter);
    expect_cholesky(Sigma);
  }

  // time-series correlation
  for (double rho : std::vector<double>{0.0, 0.9}) {
    for (size_t n = 1; n < 5; ++n) {
      Eigen::MatrixXd Sigma(n, n);
      for (int i = 0; i < n; ++i) {
        Sigma(i, i) = 1;
        for (int j = 0; j < i; ++j) {
          Sigma(i, j) = std::pow(rho, fabs(i - j));
          Sigma(j, i) = Sigma(i, j);
        }
      }
      expect_cholesky(Sigma);
    }
  }
}
