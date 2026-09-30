#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/fun/multiply_operator_util.hpp>

// Eigen products of complex autodiff matrices are evaluated by Stan's
// product hooks.

using stan::test::multiply_operator::expect_hooked;

namespace {
Eigen::MatrixXcd complex_matrix() {
  Eigen::MatrixXcd x(2, 2);
  x << std::complex<double>(1, 0.5), 2, std::complex<double>(3, -1), 4;
  return x;
}
}  // namespace

TEST(mathMixMultiplyOperator, complexMatrixProduct) {
  auto xpr = [](const auto& x, const auto& y) { return x * y; };
  auto f = [](const auto& x, const auto& y) { return (x * y).eval(); };
  Eigen::MatrixXd d_mat(2, 2);
  d_mat << 1, 2, 3, 4;
  Eigen::MatrixXcd cd_mat = complex_matrix();
  Eigen::VectorXcd cd_vec(2);
  cd_vec << std::complex<double>(0.5, 1), 2;
  expect_hooked(xpr, cd_mat, cd_mat);
  expect_hooked(xpr, cd_mat, cd_vec);
  stan::test::expect_ad(f, d_mat, cd_mat);
  stan::test::expect_ad(f, cd_mat, d_mat);
  stan::test::expect_ad(f, cd_mat, cd_mat);
  stan::test::expect_ad(f, cd_mat, cd_vec);
}

TEST(mathMixMultiplyOperator, complexCwiseAndScalar) {
  auto cwise = [](const auto& x, const auto& y) { return x.cwiseProduct(y); };
  auto scale = [](const auto& s, const auto& x) { return s * x; };
  auto f_cwise
      = [](const auto& x, const auto& y) { return x.cwiseProduct(y).eval(); };
  auto f_scale = [](const auto& s, const auto& x) { return (s * x).eval(); };
  Eigen::MatrixXd d_mat(2, 2);
  d_mat << 1, 2, 3, 4;
  Eigen::MatrixXcd cd_mat = complex_matrix();
  std::complex<double> cd_scalar(0.5, -0.25);
  expect_hooked(cwise, cd_mat, cd_mat);
  expect_hooked(scale, cd_scalar, cd_mat);
  stan::test::expect_ad(f_cwise, cd_mat, cd_mat);
  stan::test::expect_ad(f_cwise, d_mat, cd_mat);
  stan::test::expect_ad(f_scale, cd_scalar, d_mat);
  stan::test::expect_ad(f_scale, cd_scalar, cd_mat);
}
