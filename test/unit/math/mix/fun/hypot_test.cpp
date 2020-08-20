#include <test/unit/math/test_ad.hpp>

TEST(mathMixScalFun, hypot) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::hypot(x1, x2);
  };
  stan::test::expect_common_nonzero_binary(f);
  stan::test::expect_ad(f, 0.5, 0.0);
  stan::test::expect_ad(f, 0.5, 2.3);
  stan::test::expect_ad(f, 0.5001, 1.0001);
  stan::test::expect_ad(f, 3.0, 4.0);
  stan::test::expect_ad(f, 3.0, 6.0);
}

TEST(mathMixScalFun, hypot_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::hypot;
    return hypot(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 3, 1;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}
