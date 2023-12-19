#include <test/unit/math/test_ad.hpp>
#include <stdexcept>

TEST(MathMixMatFun, singularValues) {
  auto f = [](const auto& x) { return stan::math::singular_values(x); };

  Eigen::MatrixXd m00(0, 0);
  EXPECT_NO_THROW(f(m00));

  Eigen::MatrixXd m11(1, 1);
  m11 << 1.1;
  stan::test::expect_ad(f, m11);
  stan::test::expect_ad_matvar(f, m11);

  Eigen::MatrixXd m22(2, 2);
  m22 << 3, -5, 7, 11;
  stan::test::expect_ad(f, m22);
  stan::test::expect_ad_matvar(f, m22);

  Eigen::MatrixXd m23(2, 3);
  m23 << 3, 5, -7, -11, 13, -17;
  Eigen::MatrixXd m32 = m23.transpose();
  stan::test::expect_ad(f, m23);
  stan::test::expect_ad(f, m32);
  stan::test::expect_ad_matvar(f, m23);
  stan::test::expect_ad_matvar(f, m32);

  Eigen::MatrixXd a22(2, 2);
  a22 << 1, 2, 3, 4;
  stan::test::expect_ad(f, a22);
  stan::test::expect_ad_matvar(f, a22);

  Eigen::MatrixXcd c22(2, 2);
  a22 << 1, 2, 3, 4;
  stan::test::expect_ad(f, a22);
  stan::test::expect_ad_matvar(f, a22);
}
