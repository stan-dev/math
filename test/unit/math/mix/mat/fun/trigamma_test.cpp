#include <test/unit/math/test_ad.hpp>

TEST(mathMixMatFun, trigamma) {
  auto f = [](const auto& x1) { return stan::math::trigamma(x1); };
  // works, but not much precision with negative numbers
  std::vector<double> xs{0.49, stan::math::positive_infinity(),
                         stan::math::negative_infinity(),
                         stan::math::not_a_number()};
  for (double x : xs)
    stan::test::expect_ad_vectorized(f, x);
}
