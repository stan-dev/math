#include <test/unit/math/test_ad.hpp>

TEST(mixFunctor, integrate1DDoubleExponential) {
  auto f = [&](const auto& x_input, const auto& lb, const auto& ub) {
    auto func = [](const auto& x, const auto& xc, std::ostream* msgs,
                   const auto& theta) {
      return stan::math::exp(theta * stan::math::cos(2 * 3.141593 * x)) + theta;
    };
    const double relative_tolerance = std::sqrt(stan::math::EPSILON);
    const double absolute_tolerance = 0.0;
    const int max_refinements = 15;
    std::ostringstream* msgs = nullptr;
    return stan::math::integrate_1d_double_exponential_tol(
        func, lb, ub, relative_tolerance, absolute_tolerance, max_refinements,
        msgs, x_input);
  };
  stan::test::expect_ad(f, 0.75, 0, 1);
  stan::test::expect_ad(f, 0.2, 0.2, 0.7);
  // The NOT_A_NUMBER case is included in the upstream integrate_1d mix
  // test. We keep it here too because tanh_sinh's behaviour on
  // NaN-everywhere integrands matches what integrate_1d_double_exponential
  // inherits (the upstream test passes; this is a port of the same).
  stan::test::expect_ad(f, stan::math::NOT_A_NUMBER, 0, 1);
}
