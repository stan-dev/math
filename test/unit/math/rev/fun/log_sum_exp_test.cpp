#include <stan/math/prim.hpp>
#include <stan/math/rev.hpp>
#include <test/unit/math/rev/expect_identity.hpp>
#include <gtest/gtest.h>

TEST(MathFunctions, log_sum_exp_identities_rev) {
  using stan::math::log_sum_exp;
  using stan::math::var;
  using stan::test::expect_identity;

  std::vector<double> values_to_check = 
   {1e-8, 3e-5, 8e-1, 1, 8.6, 4.3e3, 9.6e6, 9e9, 1.8e15, 4.3e30};

  auto lh_duplicate = [](const var& x) { return log_sum_exp(x, x); };
  auto rh_duplicate = [](const var& x) { return stan::math::LOG_TWO + x; };

  auto lh_plus_one = [](const var& x) { return log_sum_exp(x, x + 1); };
  auto rh_plus_one = [](const var& x) { return x + log(1 + stan::math::e()); };

  for(double x : values_to_check) {
    expect_identity("Duplicate argument", lh_duplicate, rh_duplicate, x);
    expect_identity("Plus one", lh_plus_one, rh_plus_one, x);
  }

  auto lh_multiply = [](const var& x, const var& y) {
    return log_sum_exp(1e5 + x, 1e5 + y);
  };

  auto rh_multiply = [](const var& x, const var& y) { 
    return 1e5 + log_sum_exp(x, y); 
  };
  auto lh_orig = [](const var& x, const var& y) {
    return log_sum_exp(x, y);
  };

  auto rh_distribute = [](const var& x, const var& y) {
    return x + log_sum_exp(0, y - x); 
  };

  for(double x : values_to_check) {
    for(double y : values_to_check) {
      expect_identity("Multiply", lh_multiply, rh_multiply, x, y);
      expect_identity("Distribute", lh_orig, rh_distribute, x, y);
    }
  }  
  // log(exp(1e5))
}
