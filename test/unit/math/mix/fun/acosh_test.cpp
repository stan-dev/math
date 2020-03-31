#include <test/unit/math/test_ad.hpp>
#include <vector>

TEST(mathMixMatFun, acosh) {
  auto f = [](const auto& x1) {
    using stan::math::acosh;
    return acosh(x1);
  };
  for (double x : stan::test::internal::common_args())
    stan::test::expect_unary_vectorized(x);
  stan::test::expect_unary_vectorized(f, 1.5, 3.2, 5, 10, 12.9);
  // avoid pole at complex zero that can't be autodiffed
  for (double re : std::vector<double>{-0.2, 0, 0.3}) {
    for (double im : std::vector<double>{-0.3, 0.2}) {
      stan::test::expect_ad(f, std::complex<double>{re, im});
    }
  }
}
