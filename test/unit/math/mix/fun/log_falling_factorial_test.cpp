#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, logFallingFactorial) {
  auto f = [](const auto& x1, const auto& x2) {
    return stan::math::log_falling_factorial(x1, x2);
  };
  stan::test::expect_ad(f, -3.0, -2);    // throws
  stan::test::expect_ad(f, -3.0, 2);     // throws
  stan::test::expect_ad(f, 1.0, -3.0);   // throws
  stan::test::expect_ad(f, -3.0, -3.0);  // throws

  stan::test::expect_ad(f, 2.1, 4);
  stan::test::expect_ad(f, 4.0, 2);
  stan::test::expect_ad(f, 4.0, 4);
  stan::test::expect_ad(f, 5.0, 4.0);
  stan::test::expect_ad(f, 6.0, 4.0);

  stan::test::expect_ad(f, 5, 3);
  stan::test::expect_ad(f, 5, 3.0);
  stan::test::expect_ad(f, 5.0, 3);
  stan::test::expect_ad(f, 5.0, 3.0);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 2.0, nan);
  stan::test::expect_ad(f, nan, 4.0);
  stan::test::expect_ad(f, nan, nan);
}

TEST(mathMixScalFun, logFallingFactorial_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::log_falling_factorial;
    return log_falling_factorial(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 3, 1;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}
