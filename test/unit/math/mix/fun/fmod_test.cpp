#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, fmod) {
  auto f = [](const auto& x1, const auto& x2) {
    using std::fmod;
    using stan::math::fmod;
    return fmod(x1, x2);
  };
  stan::test::expect_ad(f, 2.0, 3.0);
  stan::test::expect_ad(f, 2.0, 4.0);
  stan::test::expect_ad(f, 2.7, 1.3);
  stan::test::expect_ad(f, 3.0, 2.0);
  stan::test::expect_ad(f, 3.0, 6.0);
  stan::test::expect_ad(f, 3.0, -6.0);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);

  stan::test::expect_value(f, 2.0, 2.0);
}

TEST(mathMixScalFun, fmod_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using std::fmod;
    using stan::math::fmod;
    return fmod(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 2.0, 2.0;
  Eigen::VectorXd in2(2);
  in2 << 3.0, 4.0;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}
