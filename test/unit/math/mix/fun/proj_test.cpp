#include <test/unit/math/test_ad.hpp>
#include <complex>
#include <vector>

TEST(mixScalFun, proj) {
  auto f = [](const auto& x) { return proj(x); };
  stan::test::expect_complex_common(f);
}
