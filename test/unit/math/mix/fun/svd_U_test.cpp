#include <test/unit/math/test_ad.hpp>
#include <stdexcept>

TEST(MathMixMatFun, svd_U) {
  auto f = [](const auto& x) { return stan::math::svd_U(x); };

  Eigen::MatrixXd m00(0, 0);
  EXPECT_THROW(f(m00), std::invalid_argument);

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
  stan::test::expect_ad(f, m23);
  stan::test::expect_ad_matvar(f, m23);

  Eigen::MatrixXd m32(3, 2);
  m32 << 1, 3, -5, 7, 9, -11;
  stan::test::expect_ad(f, m32);
  stan::test::expect_ad_matvar(f, m32);

  Eigen::MatrixXd a22(2, 2);
  a22 << 1, 2, 3, 4;
  stan::test::expect_ad(f, a22);
  stan::test::expect_ad_matvar(f, a22);
}
