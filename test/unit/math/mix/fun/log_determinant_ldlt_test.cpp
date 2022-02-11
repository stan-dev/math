#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, logDeterminantLdlt) {
  auto f = [](const auto& x) {
    auto x_sym = stan::math::multiply(0.5, x + x.transpose());
    auto y = stan::math::make_ldlt_factor(x_sym);
    return stan::math::log_determinant_ldlt(y);
  };

  Eigen::MatrixXd a00(0, 0);
  stan::test::expect_ad(f, a00);
  stan::test::expect_ad_matvar(f, a00);

  Eigen::MatrixXd a(2, 2);
  a << 3, 0, 0, 4;
  stan::test::expect_ad(f, a);
  stan::test::expect_ad_matvar(f, a);

  Eigen::MatrixXd c(2, 2);
  c << 1, 0, 0, 3;
  stan::test::expect_ad(f, c);
  stan::test::expect_ad_matvar(f, c);

  Eigen::MatrixXd b(2, 2);
  b << 2, 1, 1, 3;
  stan::test::expect_ad(f, b);
  stan::test::expect_ad_matvar(f, b);

  for (const auto& rho : std::vector<double>{0, 0.9}) {
    for (const auto& y : stan::test::ar_test_cov_matrices(1, 3, rho)) {
      stan::test::expect_ad(f, y);
      stan::test::expect_ad_matvar(f, y);
    }
  }
}
