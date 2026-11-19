#include <stan/math/rev.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>

// Values and d/dy far in the tails, where the result is representable.
// mu = 0, sigma = 1, lambda = 1, alpha = 1.  References: mpmath at 700
// digits, two routes for log Phi (erfc or Laplace's continued fraction, and
// the asymptotic series or the continued fraction).
namespace normal_family_range_ends_test {

template <typename F>
void expect_tail(F f, double y, double value, double dy) {
  stan::math::var y_v = y;
  stan::math::var lp = f(y_v);
  lp.grad();
  if (std::isinf(value)) {
    EXPECT_EQ(value, lp.val()) << "value at y = " << y;
  } else {
    EXPECT_NEAR(1.0, lp.val() / value, 1e-13) << "value at y = " << y;
  }
  EXPECT_NEAR(1.0, y_v.adj() / dy, 1e-13) << "d/dy at y = " << y;
  stan::math::recover_memory();
}

}  // namespace normal_family_range_ends_test

TEST_F(AgradRev, ProbNormalRangeEnds_normal_lcdf_lccdf) {
  using normal_family_range_ends_test::expect_tail;
  using stan::math::var;
  constexpr double inf = std::numeric_limits<double>::infinity();
  auto lcdf = [](const var& y) { return stan::math::normal_lcdf(y, 0.0, 1.0); };
  auto lccdf
      = [](const var& y) { return stan::math::normal_lccdf(y, 0.0, 1.0); };
  auto std_lcdf = [](const var& y) { return stan::math::std_normal_lcdf(y); };
  expect_tail(lcdf, -1e100, -5.000000000000000159e+199,
              1.0000000000000000159e+100);
  expect_tail(std_lcdf, -1e100, -5.000000000000000159e+199,
              1.0000000000000000159e+100);
  expect_tail(lccdf, 1e100, -5.000000000000000159e+199,
              -1.0000000000000000159e+100);
  // the value overflows, the gradient is still finite
  expect_tail(lcdf, -1.7e308, -inf, 1.6999999999999999388e+308);
  expect_tail(lccdf, 1.7e308, -inf, -1.6999999999999999388e+308);
}

TEST_F(AgradRev, ProbNormalRangeEnds_lognormal_lcdf_lccdf) {
  using normal_family_range_ends_test::expect_tail;
  using stan::math::var;
  auto lcdf
      = [](const var& y) { return stan::math::lognormal_lcdf(y, 0.0, 1.0); };
  auto lccdf
      = [](const var& y) { return stan::math::lognormal_lccdf(y, 0.0, 1.0); };
  expect_tail(lcdf, std::exp(-700.0), -245007.470020909054,
              7.0996388721152466927e+306);
  expect_tail(lccdf, std::exp(700.0), -245007.47002090905398,
              -6.9017876658265548454e-302);
}

TEST_F(AgradRev, ProbNormalRangeEnds_exp_mod_normal) {
  using normal_family_range_ends_test::expect_tail;
  using stan::math::var;
  auto lpdf = [](const var& y) {
    return stan::math::exp_mod_normal_lpdf(y, 0.0, 1.0, 1.0);
  };
  auto lcdf = [](const var& y) {
    return stan::math::exp_mod_normal_lcdf(y, 0.0, 1.0, 1.0);
  };
  auto lccdf = [](const var& y) {
    return stan::math::exp_mod_normal_lccdf(y, 0.0, 1.0, 1.0);
  };
  auto cdf = [](const var& y) {
    return stan::math::exp_mod_normal_cdf(y, 0.0, 1.0, 1.0);
  };
  expect_tail(lpdf, -1e100, -5.000000000000000159e+199,
              1.0000000000000000159e+100);
  expect_tail(lcdf, -1e100, -5.000000000000000159e+199,
              1.0000000000000000159e+100);
  expect_tail(lcdf, -30.0, -457.75733801210989603, 30.065381239817942414);
  expect_tail(lccdf, -30.0, -1.5794797074175364879e-199,
              -4.7487659564064334107e-198);
  expect_tail(lccdf, 1e300, -1.0000000000000000525e+300, -1.0);
  expect_tail(cdf, -30.0, 1.5794797074175364879e-199,
              4.7487659564064334107e-198);
}

TEST_F(AgradRev, ProbNormalRangeEnds_skew_normal_lpdf) {
  using normal_family_range_ends_test::expect_tail;
  using stan::math::var;
  constexpr double inf = std::numeric_limits<double>::infinity();
  auto lpdf = [](const var& y) {
    return stan::math::skew_normal_lpdf(y, 0.0, 1.0, 1.0);
  };
  expect_tail(lpdf, -1e100, -1.0000000000000000318e+200,
              2.0000000000000000318e+100);
  // the value overflows, the gradient is still finite
  expect_tail(lpdf, -1e300, -inf, 2.000000000000000105e+300);
}
