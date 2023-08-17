#include <test/unit/math/test_ad.hpp>
#include <stdexcept>

TEST(MathMixMatFun, svd) {
  auto f = [](const auto& x) { return std::get<0>(stan::math::svd(x)); };
  auto g = [](const auto& x) { return std::get<1>(stan::math::svd(x)); };
  auto h = [](const auto& x) { return std::get<2>(stan::math::svd(x)); };

  Eigen::MatrixXd m00(0, 0);
  EXPECT_NO_THROW(f(m00));
  EXPECT_NO_THROW(g(m00));
  EXPECT_NO_THROW(h(m00));

  Eigen::MatrixXd m11(1, 1);
  m11 << 1.1;
  stan::test::expect_ad(f, m11);
  stan::test::expect_ad_matvar(f, m11);
  stan::test::expect_ad(g, m11);
  stan::test::expect_ad_matvar(g, m11);
  stan::test::expect_ad(h, m11);
  stan::test::expect_ad_matvar(h, m11);

  Eigen::MatrixXd m22(2, 2);
  m22 << 3, -5, 7, 11;
  stan::test::expect_ad(f, m22);
  stan::test::expect_ad_matvar(f, m22);
  stan::test::expect_ad(g, m22);
  stan::test::expect_ad_matvar(g, m22);
  stan::test::expect_ad(h, m22);
  stan::test::expect_ad_matvar(h, m22);

  Eigen::MatrixXd m23(2, 3);
  m23 << 3, 5, -7, -11, 13, -17;
  stan::test::expect_ad(f, m23);
  stan::test::expect_ad_matvar(f, m23);
  stan::test::expect_ad(g, m23);
  stan::test::expect_ad_matvar(g, m23);
  stan::test::expect_ad(h, m23);
  stan::test::expect_ad_matvar(h, m23);

  Eigen::MatrixXd m32(3, 2);
  m32 << 1, 3, -5, 7, 9, -11;
  stan::test::expect_ad(f, m32);
  stan::test::expect_ad_matvar(f, m32);
  stan::test::expect_ad(g, m32);
  stan::test::expect_ad_matvar(g, m32);
  stan::test::expect_ad(h, m32);
  stan::test::expect_ad_matvar(h, m32);

  Eigen::MatrixXd a22(2, 2);
  a22 << 1, 2, 3, 4;
  stan::test::expect_ad(f, a22);
  stan::test::expect_ad_matvar(f, a22);
  stan::test::expect_ad(g, a22);
  stan::test::expect_ad_matvar(g, a22);
  stan::test::expect_ad(h, a22);
  stan::test::expect_ad_matvar(h, a22);

  Eigen::MatrixXcd c22(2, 2);
  a22 << 1, 2, 3, 4;
  stan::test::expect_ad(f, a22);
  stan::test::expect_ad_matvar(f, a22);
  stan::test::expect_ad(g, a22);
  stan::test::expect_ad_matvar(g, a22);
  stan::test::expect_ad(h, a22);
  stan::test::expect_ad_matvar(h, a22);
}
