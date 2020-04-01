#include <test/unit/math/test_ad.hpp>
#include <cmath>
#include <limits>
#include <vector>

template <typename T>
void expect_scalbn(double x) {
  using std::scalbn;
  T xt(x);
  for (int n = 2; n < 5; ++n) {
    stan::test::expect_near_rel("scalbn", scalbn(x, n), scalbn(xt, n));
  }
}

void expect_all_scalbn(double x) {
  using stan::math::fvar;
  using stan::math::var;
  expect_scalbn<double>(x);
  expect_scalbn<var>(x);
  expect_scalbn<fvar<double>>(x);
  expect_scalbn<fvar<fvar<double>>>(x);
  expect_scalbn<fvar<var>>(x);
  expect_scalbn<fvar<fvar<var>>>(x);
}

TEST(mathMixMatFun, scalbn) {
  double inf = std::numeric_limits<double>::infinity();
  double nan = std::numeric_limits<double>::quiet_NaN();
  for (double x : std::vector<double>{2.3}) {
    expect_all_scalbn(x);
  }
}
