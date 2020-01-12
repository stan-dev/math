#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <boost/math/special_functions/beta.hpp>
#include <stan/math/prim/scal/fun/boost_policy.hpp>

auto f = [](const auto& a, const auto& b, const auto& x) {
  return stan::math::inc_beta(a, b, x);
};

TEST(mathMixScalFun, inc_beta) {
  // replicate pre-existing test cases
  stan::test::expect_ad(f, 0.6, 0.3, 0.5);
  stan::test::expect_ad(f, 0.9, 0.9, 0.9);
  stan::test::expect_ad(f, 1.0, 1.0, 0.4);
}

TEST(mathMixScalFun, inc_beta_integer) {
  stan::test::expect_ad(f, 3, 2, 0.2);
  stan::test::expect_ad(f, 3, 2.0, 0.2);
  stan::test::expect_ad(f, 3.0, 2, 0.2);
  stan::test::expect_ad(f, 3.0, 2.0, 0.2);
}

TEST(mathMixScalFun, inc_beta_nan) {
  double nan = stan::math::NOT_A_NUMBER;

  stan::test::expect_ad(f, 0.6, 0.3, nan);
  stan::test::expect_ad(f, 0.6, nan, 0.5);
  stan::test::expect_ad(f, 0.6, nan, nan);
  stan::test::expect_ad(f, nan, 0.3, 0.5);
  stan::test::expect_ad(f, nan, 0.3, nan);
  stan::test::expect_ad(f, nan, nan, 0.5);
  stan::test::expect_ad(f, nan, nan, nan);
}

TEST(mathMixScalFun, inc_beta_a_boundary) {
  const double inf = stan::math::INFTY;
  stan::test::expect_ad(f, inf, 0.5, 0.5);
  stan::test::expect_all_throw(f, -0.1, 0.5, 0.5);
}

TEST(mathMixScalFun, inc_beta_b_boundary) {
  const double inf = stan::math::INFTY;
  stan::test::expect_ad(f, 0.5, inf, 0.5);
  stan::test::expect_all_throw(f, 0.5, -0.1, 0.5);
}

TEST(mathMixScalFun, inc_beta_x_boundary) {
  // Don't try to finite difference right next to the boundaries, the perturbation
  // goes across and causes spurious failures.
  stan::test::expect_value(f, 0.5, 0.5, 0);
  stan::test::expect_value(f, 0.5, 0.5, 1);
  stan::test::expect_ad(f, 0.5, 0.5, -0.1);
  stan::test::expect_ad(f, 0.5, 0.5, 1.1);
}

TEST(mathMixScalFun, inc_beta_a_b_boundary) {
  stan::test::expect_ad(f, 0.0, 0.0, 0.5);
}
