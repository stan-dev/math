#include <stan/math/prim/scal.hpp>
#include <stan/math/prim/scal/fun/lgamma_stirling_diff.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <test/unit/math/expect_near_rel.hpp>

TEST(MathFunctions, lgamma_stirling_diff_errors) {
  // TODO[martinmodrak] nan, negative values, ...
}

TEST(MathFunctions, lgamma_stirling_diff_accuracy) {
  using stan::math::internal::lgamma_stirling_diff_big;
  using stan::math::lgamma_stirling_diff;
  using stan::test::expect_near_rel;

  long double start = std::nextafter(10, 11);
  for (long double x = start; x < lgamma_stirling_diff_big; x *= 1.5) {
  //for (long double x = start; x < 1e150; x *= 1.5) {
  //for (double x = start; x < 100; x *= 1.5) {
    long double stirling = 0.5 
      * std::log(2 * static_cast<long double>(stan::math::pi())) 
      + (x - 0.5) * std::log(x) - x;
    long double lgamma_res = std::lgamma(x);
    double diff_actual = static_cast<double>(lgamma_res - stirling);
    double diff = lgamma_stirling_diff(x);

    std::ostringstream msg;
    msg << "x = " << x << "; lgamma = " << lgamma_res
        << "; stirling = " << stirling;
    expect_near_rel(msg.str(), diff, diff_actual, 1e-4);
  }

  double before_big = std::nextafter(lgamma_stirling_diff_big, 0);
  double after_big = std::nextafter(lgamma_stirling_diff_big,
                                    stan::math::positive_infinity());
  expect_near_rel("big cutoff", lgamma_stirling_diff(before_big),
                  lgamma_stirling_diff(after_big));
}
