#include <stan/math/prim/scal/fun/selfstart_weibull.hpp>
#include <gtest/gtest.h>
#include <vector>

std::vector<double> test_selfstart_weibull(const std::vector<double>& x,
                                          double Asym, double Drop, double lrc,
                                          double pwr) {
  std::vector<double> value(x.size());

  double exp_lrc = std::exp(lrc);

  for (std::size_t i = 0; i < x.size(); ++i) {
    value[i] = Asym - Drop * std::exp(-exp_lrc * std::pow(x[i], pwr));
  }

  return value;
}

TEST(MathFunctions, selfstart_weibull) {
  using stan::math::selfstart_weibull;

  double Asym = 160;
  double Drop = 115;
  double lrc = -5.5;
  double pwr = 2;

  std::vector<double> input;
  input.push_back(0);
  input.push_back(-1);
  input.push_back(3.14);

  std::vector<double> test_results =
    test_selfstart_weibull(input, Asym, Drop, lrc, pwr);

  std::vector<double> fun_results =
    selfstart_weibull(input, Asym, Drop, lrc, pwr);

  EXPECT_EQ(test_results.size(), fun_results.size());

  for (std::size_t i = 0; i < test_results.size(); ++i)
    EXPECT_FLOAT_EQ(test_results[i], fun_results[i]);
}
