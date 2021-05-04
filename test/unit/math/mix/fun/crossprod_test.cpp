#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, crossprod) {
  auto f = [](const auto& y) { return stan::math::crossprod(y); };

  Eigen::MatrixXd t(0, 0);

  Eigen::MatrixXd x(1, 1);
  x << 3;

  Eigen::MatrixXd u(1, 3);
  u << 1, 2, 3;

  Eigen::MatrixXd y(2, 2);
  y << 3, 0, 4, -3;

  Eigen::MatrixXd v(2, 3);
  v << 1, 2, 3, -1, 4, 9;

  Eigen::MatrixXd w(3, 2);
  w << 1, 2, 3, -1, 4, 9;

  Eigen::MatrixXd z(3, 3);
  z << 1, 0, 0, 2, 3, 0, 4, 5, 6;

  for (const auto& a : std::vector<Eigen::MatrixXd>{t, x, u, y, v, w, z})
    stan::test::expect_ad(f, a);
}
