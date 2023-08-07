#include <test/unit/math/test_ad.hpp>

TEST(mixFunctor, integrate1D) {
  auto f = [&](const auto& x_input) {
    auto func = [](const auto& x, const auto& xc, std::ostream* msgs,
                    const auto& theta) {
      return stan::math::exp(theta * stan::math::cos(2 * 3.141593 * x)) + theta;
    };
    const double relative_tolerance = std::sqrt(stan::math::EPSILON);
    std::ostringstream *msgs = nullptr;
    return stan::math::integrate_1d_impl(func, 0, 1, relative_tolerance, msgs,
                                          x_input);
  };
  stan::test::expect_ad(f, 0.75);
  stan::test::expect_ad(f, 0.2);
  stan::test::expect_ad(f, stan::math::NOT_A_NUMBER);
}
