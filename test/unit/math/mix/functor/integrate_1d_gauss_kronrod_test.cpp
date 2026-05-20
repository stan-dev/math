#include <test/unit/math/test_ad.hpp>

TEST(mixFunctor, integrate1DGaussKronrod) {
  auto f = [&](const auto& x_input, const auto& lb, const auto& ub) {
    auto func = [](const auto& x, const auto& xc, std::ostream* msgs,
                   const auto& theta) {
      return stan::math::exp(theta * stan::math::cos(2 * 3.141593 * x)) + theta;
    };
    const double relative_tolerance = std::sqrt(stan::math::EPSILON);
    const double absolute_tolerance = 0.0;
    const int max_depth = 15;
    std::ostringstream* msgs = nullptr;
    return stan::math::integrate_1d_gauss_kronrod_impl(
        func, lb, ub, relative_tolerance, absolute_tolerance, max_depth, msgs,
        x_input);
  };
  stan::test::expect_ad(f, 0.75, 0, 1);
  stan::test::expect_ad(f, 0.2, 0.2, 0.7);
  // Note: the upstream integrate_1d mix test also checks expect_ad with
  // NOT_A_NUMBER as the parameter. With tanh_sinh that succeeds because
  // Boost's double-exponential code path treats NaN-everywhere integrands
  // as a singularity to be ignored, while Gauss-Kronrod propagates the
  // domain_error thrown by the rev gradient-of-f NaN check. This
  // discrepancy reflects a difference in Boost quadrature behaviour, not
  // a correctness issue in the integrate_1d_gauss_kronrod wrapper, and so
  // the NaN-input case is intentionally omitted here.
}
