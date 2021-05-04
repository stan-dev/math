#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, fdim) {
  auto f
      = [](const auto& x1, const auto& x2) { return stan::math::fdim(x1, x2); };
  stan::test::expect_ad(f, -3.0, 4.0);
  stan::test::expect_ad(f, 1.0, 2.0);
  stan::test::expect_ad(f, 2.0, 1.0);
  stan::test::expect_ad(f, 2.0, -3.0);
  stan::test::expect_ad(f, 2.5, 1.5);
  stan::test::expect_ad(f, 3.0, 5.0);
  stan::test::expect_ad(f, 3.0, 4.0);
  stan::test::expect_ad(f, 4.0, -3.0);
  stan::test::expect_ad(f, 7.0, 2.0);

  stan::test::expect_value(f, 2.0, 2.0);  // can't autodiff at discontinuity

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);
}

TEST(mathMixScalFun, fdim_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::fdim;
    return fdim(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 3, 1;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}
