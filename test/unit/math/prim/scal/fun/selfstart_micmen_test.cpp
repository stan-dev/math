#include <stan/math/prim/scal/fun/selfstart_micmen.hpp>
#include <gtest/gtest.h>
#include <vector>

std::vector<double> test_selfstart_micmen(const std::vector<double>& x,
                                          double Vm, double K) {
  std::vector<double> value(x.size());

  for (std::size_t i = 0; i < x.size(); ++i) {
    value[i] = Vm * x[i]/(K + x[i]);
  }
  return value;
}

TEST(MathFunctions, selfstart_micmen) {
  using stan::math::selfstart_micmen;

  double Vm = 200;
  double K = 0.05;

  std::vector<double> input;
  input.push_back(0);
  input.push_back(-1);
  input.push_back(3.14);

  std::vector<double> test_results = test_selfstart_micmen(input, Vm, K);

  std::vector<double> fun_results = selfstart_micmen(input, Vm, K);

  EXPECT_EQ(test_results.size(), fun_results.size());

  for (std::size_t i = 0; i < test_results.size(); ++i)
    EXPECT_FLOAT_EQ(test_results[i], fun_results[i]);
}
