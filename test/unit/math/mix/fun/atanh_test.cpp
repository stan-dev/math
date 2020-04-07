#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(mathMixMatFun, atanh) {
  auto f = [](const auto& x1) {
    using stan::math::atanh;
    return atanh(x1);
  };
  stan::test::expect_common_unary_vectorized(f);
  stan::test::expect_unary_vectorized(f, -0.9, 0.5);
  for (double re : std::vector<double>{-0.2, 0, 0.3}) {
    for (double im : std::vector<double>{-0.3, 0, 0.2}) {
      stan::test::expect_ad(f, std::complex<double>{re, im});
    }
  }
}
