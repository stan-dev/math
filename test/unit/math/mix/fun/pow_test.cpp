#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <limits>

template <typename T>
void expect_arith_instantiate() {
  using stan::math::pow;
  using std::pow;
  auto a = pow(T(1.0), 1);
  auto b = pow(T(1.0), 1.0);
  auto c = pow(1, T(1.0));
  auto d = pow(1.0, T(1.0));
}

// this one's been tricky to instantiate, so test all instances
TEST(mathMixScalFun, powInstantiations) {
  using stan::math::fvar;
  using stan::math::var;
  expect_arith_instantiate<int>();
  expect_arith_instantiate<double>();
  expect_arith_instantiate<var>();
  expect_arith_instantiate<fvar<double>>();
  expect_arith_instantiate<fvar<fvar<double>>>();
  expect_arith_instantiate<fvar<var>>();
  expect_arith_instantiate<fvar<fvar<var>>>();
}

TEST(mathMixScalFun, pow) {
  auto f = [](const auto& x1, const auto& x2) {
    using std::pow;
    using stan::math::pow;
    return pow(x1, x2);
  };
  stan::test::expect_ad(f, -0.4, 0.5);
  stan::test::expect_ad(f, 0.5, 0.5);
  stan::test::expect_ad(f, 0.5, 1.0);
  stan::test::expect_ad(f, 0.5, 1.2);
  stan::test::expect_ad(f, 0.5, 5.0);
  stan::test::expect_ad(f, 1.0, 2.0);
  stan::test::expect_ad(f, 3.0, 4.0);
  for (double y = -2; y <= 4; y += 0.5)
    stan::test::expect_ad(f, 4.0, y);

  double nan = std::numeric_limits<double>::quiet_NaN();
  stan::test::expect_ad(f, 1.0, nan);
  stan::test::expect_ad(f, nan, 1.0);
  stan::test::expect_ad(f, nan, nan);
}
