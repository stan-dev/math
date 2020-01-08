#include <test/unit/math/test_ad.hpp>
#include <limits>

template <typename T, typename U>
stan::return_type_t<T, U> g1(const T& x, const U& lb) {
  return stan::math::lb_constrain(x, lb);
}
template <typename T, typename U>
stan::return_type_t<T, U> g2(const T& x, const U& lb) {
  T lp = 0;
  return stan::math::lb_constrain(x, lb, lp);
}

template <typename T, typename U>
stan::return_type_t<T, U> g3(const T& x, const U& lb) {
  T lp = 0;
  stan::math::lb_constrain(x, lb, lp);
  return lp;
}

void expect_lb_constrain(double x, double lb) {
  auto f1 = [](const auto& x, const auto& lb) { return g1(x, lb); };
  auto f2 = [](const auto& x, const auto& lb) { return g2(x, lb); };
  auto f3 = [](const auto& x, const auto& lb) { return g3(x, lb); };
  stan::test::expect_ad(f1, x, lb);
  stan::test::expect_ad(f2, x, lb);
  stan::test::expect_ad(f3, x, lb);
}

TEST(mathMixScalFun, lbConstrain) {
  expect_lb_constrain(-1, 2);
  expect_lb_constrain(2, 4);
}
