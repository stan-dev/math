#include <stan/math/prim/scal.hpp>
#include <stan/math/prim/scal/fun/lgamma_stirling_diff.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <test/unit/math/expect_near_rel.hpp>

TEST(MathFunctions, lgamma_stirling_diff_errors) {
    //TODO[martinmodrak] nan, negative values, ...
}

TEST(MathFunctions, lgamma_stirling_diff_accuracy) {
  using stan::math::lgamma_stirling_diff;
  using stan::test::expect_near_rel;

  double start = std::nextafter(10, 11);
  for(double x = start; x < 1e150; x *= 1.5) {
      double stirling = x * (log(x) - 1) + log(x) 
        + 0.5 * (stan::math::LOG_SQRT_PI + stan::math::LOG_TWO - log(x));
      double lgamma_res = stan::math::lgamma(x);
      double diff = lgamma_stirling_diff(x);

      std::ostringstream msg;
      msg << "x = " << x;
      expect_near_rel(msg.str(), stirling + diff, lgamma_res);
  }
}

