#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, logDeterminant) {
  auto f = [](const auto& x) { return stan::math::log_determinant(x); };

  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, m00);
  stan::test::expect_ad_matvar(f, m00);

  // repeat symmetric pos def tests
  Eigen::MatrixXd a(2, 2);
  a << 3, 0, 0, 4;
  stan::test::expect_ad(f, a);
  stan::test::expect_ad_matvar(f, a);

  Eigen::MatrixXd b(2, 2);
  b << 2, 1, 1, 3;
  stan::test::expect_ad(f, b);
  stan::test::expect_ad_matvar(f, b);

  Eigen::MatrixXd c(2, 2);
  c << 1, 0, 0, 3;
  stan::test::expect_ad(f, c);
  stan::test::expect_ad_matvar(f, c);

  for (const auto& rho : std::vector<double>{0, 0.9}) {
    for (const auto& y : stan::test::ar_test_cov_matrices(1, 3, rho)) {
      stan::test::expect_ad(f, y);
      stan::test::expect_ad_matvar(f, y);
    }
  }

  // asymmetric tests
  Eigen::MatrixXd d(2, 3);
  d << 1, 2, 3, 4, 5, 6;
  stan::test::expect_ad(f, d);
  stan::test::expect_ad_matvar(f, d);

  Eigen::MatrixXd e(2, 2);
  e << 0, 1, 2, 3;
  stan::test::expect_ad(f, e);
  stan::test::expect_ad_matvar(f, e);

  Eigen::MatrixXd g(2, 2);
  g << 2, 3, 5, 7;
  stan::test::expect_ad(f, g);
  stan::test::expect_ad_matvar(f, g);
}
