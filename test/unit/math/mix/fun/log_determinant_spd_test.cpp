#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, logDeterminantSpd) {
  auto f = [](const auto& x) {
    auto z = ((x + x.transpose()) * 0.5).eval();
    return stan::math::log_determinant_spd(z);
  };

  // for testing error conditions
  auto g = [](const auto& x) { return stan::math::log_determinant_spd(x); };

  Eigen::MatrixXd m00(0, 0);
  stan::test::expect_ad(f, m00);

  Eigen::MatrixXd a(2, 2);
  a << 3, 0, 0, 4;
  stan::test::expect_ad(f, a);

  Eigen::MatrixXd b(2, 2);
  b << 2, 1, 1, 3;
  stan::test::expect_ad(f, b);

  Eigen::MatrixXd c(2, 2);
  c << 1, 0, 0, 3;
  stan::test::expect_ad(f, c);

  for (const auto& rho : std::vector<double>{0, 0.9}) {
    for (const auto& y : stan::test::ar_test_cov_matrices(1, 3, rho)) {
      stan::test::expect_ad(f, y);
    }
  }

  Eigen::MatrixXd d(2, 3);  // not square
  d << 1, 2, 3, 4, 5, 6;
  stan::test::expect_ad(g, d);

  Eigen::MatrixXd e(2, 2);  // asymmetric
  e << 1, 2, 3, 4;
  stan::test::expect_ad(g, e);
}
