#include <test/unit/math/test_ad.hpp>
#include <limits>

TEST(mathMixScalFun, fmax) {
  auto f
      = [](const auto& x1, const auto& x2) { return stan::math::fmax(x1, x2); };
  stan::test::expect_ad(f, -3.0, 4.0);
  stan::test::expect_ad(f, 1.3, 2.0);
  stan::test::expect_ad(f, 2.3, 2.0);
  stan::test::expect_ad(f, 2.5, 1.5);
  stan::test::expect_ad(f, 4.0, -3.0);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);

  stan::test::expect_value(f, 2.0, 2.0);
}

TEST(mathMixScalFun, fmax_vec) {
  auto f = [](const auto& x1, const auto& x2) {
    using stan::math::fmax;
    return fmax(x1, x2);
  };

  Eigen::VectorXd in1(2);
  in1 << 3, 1;
  Eigen::VectorXd in2(2);
  in2 << 0.5, 3.4;
  stan::test::expect_ad_vectorized_binary(f, in1, in2);
}

TEST(mathMixScalFun, fmax_equal) {
  using stan::math::fmax;
  using stan::math::fvar;
  using stan::math::var;

  var a = 1, b = 1, c = fmax(a, b);
  c.grad();

  fvar<double> d = {1, 1}, e = {1, 1}, f = fmax(d, e);
  EXPECT_FLOAT_EQ(a.adj() + b.adj(), f.d_);
}
