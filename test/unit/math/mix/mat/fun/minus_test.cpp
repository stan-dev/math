#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(MathMixMatFun, minus) {
  auto f = [](const auto& x) { return stan::math::minus(x); };

  double x = 10;
  stan::test::expect_ad(f, x);

  Eigen::VectorXd y0(0);
  Eigen::VectorXd y1(1);
  y1 << -100;
  Eigen::VectorXd y2(2);
  y2 << -100, 0;
  Eigen::VectorXd y3(3);
  y3 << -100, 0, 1;
  for (const auto& y : std::vector<Eigen::VectorXd>{y0, y1, y2, y3}) {
    stan::test::expect_ad(f, y);
  }

  Eigen::RowVectorXd z0(0);
  Eigen::RowVectorXd z1(1);
  z1 << -100;
  Eigen::RowVectorXd z2(2);
  z2 << -100, 0;
  Eigen::RowVectorXd z3(3);
  z3 << -100, 0, 1;
  for (const auto& y : std::vector<Eigen::RowVectorXd>{z0, z1, z2, z3}) {
    stan::test::expect_ad(f, y);
  }

  Eigen::MatrixXd u00(0, 0);
  Eigen::MatrixXd u11(1, 1);
  u11 << 1;
  Eigen::MatrixXd u22(2, 2);
  u22 << 1, 2, 3, 4;
  Eigen::MatrixXd u23(2, 3);
  u23 << -100, 0, 1, 20, -40, 2;
  for (const auto& y : std::vector<Eigen::MatrixXd>{u00, u11, u22, u23}) {
    stan::test::expect_ad(f, y);
  }
}
