#include <test/unit/math/test_ad.hpp>
#include <limits>

template <typename T>
auto g1(const T& x) {
  auto x_cons = stan::math::identity_constrain(x);
  auto x_free = stan::math::identity_free(x_cons);
  return x_free;
}
template <typename T>
auto g2(const T& x) {
  stan::value_type_t<T> lp = 0;
  auto x_cons = stan::math::identity_constrain(x, lp);
  auto x_free = stan::math::identity_free(x_cons);
  return x_free;
}
template <typename T>
auto g3(const T& x) {
  stan::value_type_t<T> lp = 0;
  auto x_cons = stan::math::identity_constrain(x, lp);
  auto x_free = stan::math::identity_free(x_cons);
  return lp;
}
template <typename T>
void expect_identity_constrain(const T& x) {
  auto f1 = [](const auto& x) { return g1(x); };
  auto f2 = [](const auto& x) { return g2(x); };
  auto f3 = [](const auto& x) { return g3(x); };
  stan::test::expect_ad(f1, x);
  stan::test::expect_ad(f2, x);
  stan::test::expect_ad(f3, x);
}

TEST(mathMixScalFun, identityConstrain) {
  expect_identity_constrain(-1);
  expect_identity_constrain(2.1);
}
