#include <stan/math/prim.hpp>
#include <stan/math/rev.hpp>
#include <test/unit/math/rev/expect_identity.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, log_sum_exp_identities_rev) {
  using stan::test::expect_identity;
  using stan::math::log_sum_exp;
  using stan::math::var;

  auto lh0 = []() {
    return var(3);
  };

  auto rh0 = []() {
    return var(3);
  };

  expect_identity(lh0, rh0);

  auto lh1 = [](const var& x) {
    return log_sum_exp(x, x);
  };

  auto rh1 = [](const var& x) {
    return stan::math::LOG_TWO + x;
  };

  expect_identity(lh1, rh1, 1);

  auto lh2 = [](const var& x, const var& y) {
    return log_sum_exp(1e5 + x, 1e5 + y);
  };

  auto rh2 = [](const var& x, const var& y) {
    return 1e5 + log_sum_exp(x, y);
  };

  expect_identity(lh2, rh2, 1, 2);
  // log(exp(1e5))
}
