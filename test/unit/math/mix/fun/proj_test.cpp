#include <test/unit/math/test_ad.hpp>
#include <test/unit/math/mix/util.hpp>
#include <complex>
#include <vector>

TEST_F(mathMix,  proj) {
  auto f = [](const auto& x) { return proj(x); };
  stan::test::expect_complex_common(f);
}
