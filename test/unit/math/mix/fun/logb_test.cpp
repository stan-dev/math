#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <limits>
#include <vector>

template <typename T>
void expect_logb(double x) {
  using std::logb;
  T xt(x);
  stan::test::expect_near_rel("logb", logb(x), logb(xt));
}

void expect_all_logb(double x) {
  using stan::math::fvar;
  using stan::math::var;
  expect_logb<double>(x);
  expect_logb<var>(x);
  expect_logb<fvar<double>>(x);
  expect_logb<fvar<fvar<double>>>(x);
  expect_logb<fvar<var>>(x);
  expect_logb<fvar<fvar<var>>>(x);
}

TEST(mathMixMatFun, logb) {
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  for (double x : std::vector<double>{-inf, -2.3, 0, 1, 2, 3.9, inf, nan}) {
    expect_all_logb(x);
  }
}
