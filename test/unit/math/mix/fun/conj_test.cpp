#include <test/unit/math/test_ad.hpp>
#include <complex>
#include <vector>

TEST(mixScalFun, conj) {
  auto f = [](const auto& x) { return conj(x); };
  stan::test::expect_complex_common(f);
}
