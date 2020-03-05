#include <test/unit/math/test_ad.hpp>
#include <complex>
#include <vector>

TEST(mixScalFun, arg) {
  auto f = [](const auto& x) { return norm(x); };
  stan::test::expect_complex_common(f);
}
