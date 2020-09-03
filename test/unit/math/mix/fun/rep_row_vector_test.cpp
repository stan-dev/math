#include <test/unit/math/test_ad.hpp>

TEST(MathMixMatFun, repRowVector) {
  auto f = [](int n) {
    return [=](const auto& y) { return stan::math::rep_row_vector(y, n); };
  };
  double y = 3;
  stan::test::expect_ad(f(0), y);
  stan::test::expect_ad(f(1), y);
  stan::test::expect_ad(f(4), y);
}
