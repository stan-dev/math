#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <vector>
/*
TEST(MathMixMatFun, determinant) {
  using stan::test::relative_tolerance;
  auto f = [](const auto& y) { return stan::math::determinant(y); };

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = relative_tolerance(1e-4, 1e-2);
  tols.hessian_fvar_hessian_ = relative_tolerance(1e-4, 1e-2);

  Eigen::MatrixXd z(0, 0);

  Eigen::MatrixXd a(1, 1);
  a << -1;

  Eigen::MatrixXd b(2, 2);
  b << 2, 3, 5, 7;

  Eigen::MatrixXd c(2, 2);
  c << 1, 0.9, 0.9, 1;

  Eigen::MatrixXd d(3, 3);
  d << 1, 2, 3, 13, 17, 19, 23, 11, 7;

  Eigen::MatrixXd e(2, 2);
  e << 0, 1, 2, 3;

  Eigen::MatrixXd g(4, 4);
  for (int i = 0; i < 4; ++i) {
    g(i, i) = 1;
    for (int j = 0; j < i; ++j) {
      g(i, j) = std::pow(0.9, std::fabs(i - j));
      g(j, i) = g(i, j);
    }
  }

  for (const auto& x : std::vector<Eigen::MatrixXd>{z, a, b, c, d, e, g}) {
    stan::test::expect_ad(tols, f, x);
    stan::test::expect_ad_matvar(f, x);
  }
}
*/
TEST(MathMixMatFun, determinant) {
  using stan::test::relative_tolerance;
  auto f = [](const auto& y) { return stan::math::determinant(y); };

  stan::test::ad_tolerances tols;
  tols.hessian_hessian_ = relative_tolerance(1e-4, 1e-2);
  tols.hessian_fvar_hessian_ = relative_tolerance(1e-4, 1e-2);
  using complex_d = std::complex<double>;
  using complex_matrix = Eigen::Matrix<complex_d, -1, -1>;
  complex_matrix z(0, 0);

  complex_matrix a(1, 1);
  a << complex_d(-1);

  complex_matrix b(2, 2);
  b << complex_d(2), complex_d(3), complex_d(5), complex_d(7);

  complex_matrix c(2, 2);
  c << complex_d(1), complex_d(0.9), complex_d(0.9), complex_d(1);

  complex_matrix d(3, 3);
  d << complex_d(1), complex_d(2), complex_d(3), complex_d(13), complex_d(17), complex_d(19), complex_d(23), complex_d(11), complex_d(7);

  complex_matrix e(2, 2);
  e << complex_d(0), complex_d(1), complex_d(2), complex_d(3);

  complex_matrix g(4, 4);
  for (int i = 0; i < 4; ++i) {
    g(i, i) = complex_d(1);
    for (int j = 0; j < i; ++j) {
      g(i, j) = complex_d(std::pow(0.9, std::fabs(i - j)));
      g(j, i) = g(i, j);
    }
  }

  for (const auto& x : std::vector<complex_matrix>{e}) {
    stan::test::expect_ad(tols, f, x);
  }
}
