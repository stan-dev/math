#include <test/unit/math/test_ad.hpp>
#include <limits>

namespace identity_constrain_test {
template <typename T>
T g1(const T& x) {
  return stan::math::identity_constrain(x);
}
template <typename T>
T g2(const T& x) {
  T lp = 0;
  return stan::math::identity_constrain(x, lp);
}
template <typename T>
T g3(const T& x) {
  T lp = 0;
  stan::math::identity_constrain(x, lp);
  return lp;
}
}  // namespace identity_constrain_test

void expect_identity_constrain(double x) {
  auto f1 = [](const auto& x) { return identity_constrain_test::g1(x); };
  auto f2 = [](const auto& x) { return identity_constrain_test::g2(x); };
  auto f3 = [](const auto& x) { return identity_constrain_test::g3(x); };
  stan::test::expect_ad(f1, x);
  stan::test::expect_ad(f2, x);
  stan::test::expect_ad(f3, x);
}

TEST(mathMixScalFun, identityConstrain) {
  expect_identity_constrain(-1);
  expect_identity_constrain(2.1);
}
