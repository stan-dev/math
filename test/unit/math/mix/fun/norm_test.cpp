#include <test/unit/math/test_ad.hpp>

TEST(mixScalFun, norm) {
  auto f = [](const auto& x) { return norm(x); };
  stan::test::expect_complex_common(f);
}
