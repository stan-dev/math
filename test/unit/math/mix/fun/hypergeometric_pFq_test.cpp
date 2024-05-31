#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, hyper_2f2) {
  auto f = [](const auto& a, const auto& b, const auto& z) {
    using stan::math::hypergeometric_pFq;
    return hypergeometric_pFq(a, b, z);
  };

  Eigen::VectorXd in1(2);
  in1 << 4, 2;
  Eigen::VectorXd in2(2);
  in2 << 6, 3;
  double z = 4;

  stan::test::expect_ad(f, in1, in2, z);
}

TEST(mathMixScalFun, hyper_2f3) {
  auto f = [](const auto& a, const auto& b, const auto& z) {
    using stan::math::hypergeometric_pFq;
    return hypergeometric_pFq(a, b, z);
  };

  Eigen::VectorXd in1(2);
  in1 << 2, 3;
  Eigen::VectorXd in2(3);
  in2 << 2, 4, 5;
  double z = 1;

  stan::test::expect_ad(f, in1, in2, z);
}

TEST(mathMixScalFun, hyper_4f3) {
  auto f = [](const auto& a, const auto& b, const auto& z) {
    using stan::math::hypergeometric_pFq;
    return hypergeometric_pFq(a, b, z);
  };

  Eigen::VectorXd in1(4);
  in1 << 1, 2, 3, 4;
  Eigen::VectorXd in2(3);
  in2 << 5, 6, 7;
  double z = 0.8;

  stan::test::expect_ad(f, in1, in2, z);
}

TEST(mathMixScalFun, hyper_2f1) {
  auto f = [](const auto& a, const auto& b, const auto& z) {
    using stan::math::hypergeometric_pFq;
    return hypergeometric_pFq(a, b, z);
  };

  Eigen::VectorXd in1(2);
  in1 << 1, 1;
  Eigen::VectorXd in2(1);
  in2 << 1;
  double z = 0.6;

  stan::test::expect_ad(f, in1, in2, z);

  in1 << 1, 31;
  in2 << 41;
  z = 0.6;
  stan::test::expect_ad(f, in1, in2, z);

  in1 << 1, -2.1;
  in2 << 41;
  z = 0.6;
  stan::test::expect_ad(f, in1, in2, z);

  in1 << 1, -0.5;
  in2 << 10;
  z = 0.3;
  stan::test::expect_ad(f, in1, in2, z);

  in1 << 1, -0.5;
  in2 << 10.6;
  z = 0.3;
  stan::test::expect_ad(f, in1, in2, z);

  in1 << -0.5, -4.5;
  in2 << 11;
  z = 0.3;
  stan::test::expect_ad(f, in1, in2, z);

  in1 << -0.5, -4.5;
  in2 << -3.2;
  z = 0.9;
  stan::test::expect_ad(f, in1, in2, z);

  in1 << 2, 1;
  in2 << 2;
  z = 0.4;
  stan::test::expect_ad(f, in1, in2, z);

  in1 << 3.70975, 1.0;
  in2 << 2.70975;
  z = -0.2;
  stan::test::expect_ad(f, in1, in2, z);
}

TEST(mathMixScalFun, hyper_3f2) {
  using stan::math::var;

  auto f = [](const auto& a, const auto& b, const auto& z) {
    using stan::math::hypergeometric_pFq;
    return hypergeometric_pFq(a, b, z);
  };

  Eigen::VectorXd in1(3);
  in1 << 1.0, 1.0, 1.0;
  Eigen::VectorXd in2(2);
  in2 << 1.0, 1.0;
  double z = 0.6;

  stan::test::expect_ad(f, in1, in2, z);

  in1 << 1.0, -0.5, -2.5;
  in2 << 10.0, 1.0;
  z = 0.3;
  stan::test::expect_ad(f, in1, in2, z);
}
