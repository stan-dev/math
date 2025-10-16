#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <vector>

namespace cholesky_decompose_test {
// can't autodiff directly through Cholesky due to symmetry test;
// use unconstrained input and constrain to test Cholesky derivs;
// dof must be (n choose 2) + n
auto f(int dof) {
  return [=](const auto& x) {
    stan::math::promote_scalar_t<stan::value_type_t<decltype(x)>,
                                 Eigen::Matrix<double, -1, -1>>
        y;
    try {
      y = stan::math::cov_matrix_constrain(x, dof);
    } catch (...) {
      ADD_FAILURE() << "FAILED AT COV_MATRIX_CONSTRAIN";
      throw;
    }
    return stan::math::cholesky_decompose(y);
  };
}

auto f_matvar = [](const auto& x) { return stan::math::cholesky_decompose(x); };

void expect_cholesky(const Eigen::MatrixXd& Sigma) {
  Eigen::VectorXd yy = stan::math::cov_matrix_free(Sigma);
  // lazy, solving for x in x = (N * (N + 1)) / 2
  int dof = .5 * (std::sqrt(8 * yy.size() + 1) - 1);
  stan::test::expect_ad(f(dof), yy);
}
}  // namespace cholesky_decompose_test

TEST(MathMixMatFun, choleskyDecomposeSpecific) {
  // 1 x 1 matrix;  (1 choose 2) + 1 = 1
  Eigen::VectorXd x1(1);
  x1 << 1;
  stan::test::expect_ad(cholesky_decompose_test::f(1), x1);
  Eigen::MatrixXd x1_mat = stan::math::cov_matrix_constrain(x1, 1);
  stan::test::expect_ad_matvar(cholesky_decompose_test::f_matvar, x1_mat);

  // 2 x 2 matrix;  (2 choose 2) + 2 = 3
  Eigen::VectorXd x3(3);
  x3 << 1, 2, -1;
  stan::test::expect_ad(cholesky_decompose_test::f(2), x3);
  Eigen::MatrixXd x3_mat = stan::math::cov_matrix_constrain(x3, 2);
  stan::test::expect_ad_matvar(cholesky_decompose_test::f_matvar, x3_mat);

  // 3 x 3 matrix;  (3 choose 2) + 3 = 6
  Eigen::VectorXd x6(6);
  x6 << 1, -1, 1.1, 1.4, 2.1, 0.7;
  stan::test::expect_ad(cholesky_decompose_test::f(3), x6);
  Eigen::MatrixXd x6_mat = stan::math::cov_matrix_constrain(x6, 3);
  stan::test::expect_ad_matvar(cholesky_decompose_test::f_matvar, x6_mat);

  // 4 x 4 matrix;  (4 choose 2) + 4 = 10
  Eigen::VectorXd x10(10);
  x10 << 1, -0.1, 1.1, 1.4, -1.1, 0.7, 1.0, 1.3, -0.5, 0.3;
  stan::test::expect_ad(cholesky_decompose_test::f(4), x10);
  Eigen::MatrixXd x10_mat = stan::math::cov_matrix_constrain(x10, 4);
  stan::test::expect_ad_matvar(cholesky_decompose_test::f_matvar, x10_mat);

  // 2 x 3 matrix will throw; test directly
  auto g = [](const auto& x) { return stan::math::cholesky_decompose(x); };
  Eigen::MatrixXd y(2, 3);
  y << 1, 2, 3, 4, 5, 6;
  stan::test::expect_ad(g, y);
  stan::test::expect_ad_matvar(g, y);

  // asymmetric will throw
  Eigen::MatrixXd z(2, 2);
  z << 1, 2, 3, 4;
  stan::test::expect_ad(g, z);
  stan::test::expect_ad_matvar(g, y);
}

TEST(MathMixMatFun, choleskyDecomposeGeneral) {
  // general sizes
  for (int n = 0; n < 8; ++n) {
    int dof = (n * (n + 1)) / 2;
    Eigen::VectorXd y(dof);
    for (int i = 0; i < dof; ++i)
      y(i) = (i * 10) / 10000.0;
    stan::test::ad_tolerances tol;
    using stan::test::relative_tolerance;
    stan::test::expect_ad(tol, cholesky_decompose_test::f(n), y);

    Eigen::MatrixXd y_mat = stan::math::cov_matrix_constrain(y, n);
    stan::test::expect_ad_matvar(cholesky_decompose_test::f_matvar, y_mat);
    stan::math::recover_memory();
  }
}

TEST(MathMixMatFun, choleskyDecomposeGeneralBig) {
  // general sizes
  for (int n : std::vector<int>{12, 13, 20, 40, 50}) {
    int dof = (n * (n + 1)) / 2;
    Eigen::VectorXd y(dof);
    for (int i = 0; i < dof; ++i)
      y(i) = (i * 10) / 10000.0;
    stan::test::ad_tolerances tol;
    using stan::test::relative_tolerance;
    Eigen::MatrixXd y_mat = stan::math::cov_matrix_constrain(y, n);
    stan::test::expect_ad_matvar(cholesky_decompose_test::f_matvar, y_mat);
    stan::math::recover_memory();
  }
}

// GP covar
TEST(MathMixMatFun, choleskyDecomposeGP) {
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
    cholesky_decompose_test::expect_cholesky(Sigma);
    stan::test::expect_ad_matvar(cholesky_decompose_test::f_matvar, Sigma);
  }

  // time-series correlation
  for (double rho : std::vector<double>{0.0, 0.9}) {
    for (size_t n = 1; n < 5; ++n) {
      Eigen::MatrixXd Sigma(n, n);
      for (int i = 0; i < n; ++i) {
        Sigma(i, i) = 1;
        for (int j = 0; j < i; ++j) {
          Sigma(i, j) = std::pow(
              rho, fabs(static_cast<double>(i) - static_cast<double>(j)));
          Sigma(j, i) = Sigma(i, j);
        }
      }
      cholesky_decompose_test::expect_cholesky(Sigma);
      stan::test::expect_ad_matvar(cholesky_decompose_test::f_matvar, Sigma);
    }
  }
}
