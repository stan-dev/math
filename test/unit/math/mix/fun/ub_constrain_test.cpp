#include <test/unit/math/test_ad.hpp>
#include <limits>

namespace ub_constrain_test {
template <typename T, typename U>
stan::return_type_t<T, U> g1(const T& x, const U& ub) {
  return stan::math::ub_constrain(x, ub);
}
template <typename T, typename U>
stan::return_type_t<T, U> g2(const T& x, const U& ub) {
  T lp = 0;
  return stan::math::ub_constrain(x, ub, lp);
}
template <typename T, typename U>
stan::return_type_t<T, U> g3(const T& x, const U& ub) {
  T lp = 0;
  stan::math::ub_constrain(x, ub, lp);
  return lp;
}

void expect_ub_constrain(double x, double ub) {
  auto f1 = [](const auto& x, const auto& ub) { return g1(x, ub); };
  auto f2 = [](const auto& x, const auto& ub) { return g2(x, ub); };
  auto f3 = [](const auto& x, const auto& ub) { return g3(x, ub); };
  stan::test::expect_ad(f1, x, ub);
  stan::test::expect_ad(f2, x, ub);
  stan::test::expect_ad(f3, x, ub);
}
}  // namespace ub_constrain_test

TEST(mathMixScalFun, ub_constrain) {
  ub_constrain_test::expect_ub_constrain(-1, 2);
  ub_constrain_test::expect_ub_constrain(2, 4);
}
