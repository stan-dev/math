#include <test/unit/math/test_ad.hpp>

auto f(int n) {
  return [=](const auto& y) { return stan::math::rep_vector(y, n); };
}

TEST(MathMixMatFun, repVector) {
  double y = 3;
  stan::test::expect_ad(f(0), y);
  stan::test::expect_ad(f(1), y);
  stan::test::expect_ad(f(4), y);
}
