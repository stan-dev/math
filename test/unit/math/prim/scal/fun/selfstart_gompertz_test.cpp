#include <stan/math/prim/scal/fun/selfstart_gompertz.hpp>
#include <gtest/gtest.h>
#include <vector>

std::vector<double> test_selfstart_gompertz(
                                         const std::vector<double>& x,
                                         double Asym, double b2, double b3) {
  std::vector<double> value(x.size());

  for (std::size_t i = 0; i < x.size(); ++i) {
    value[i] = Asym * std::exp(-b2 * std::pow(b3, x[i]));
  }
  return value;
}

TEST(MathFunctions, selfstart_gompertz) {
  using stan::math::selfstart_gompertz;

  double Asym = 4.5;
  double b2 = 2.3;
  double b3 = 0.7;

  std::vector<double> input;
  input.push_back(0);
  input.push_back(-1);
  input.push_back(3.14);

  std::vector<double> test_results =
    test_selfstart_gompertz(input, Asym, b2, b3);

  std::vector<double> fun_results =
    selfstart_gompertz(input, Asym, b2, b3);

  EXPECT_EQ(test_results.size(), fun_results.size());

  for (std::size_t i = 0; i < test_results.size(); ++i)
    EXPECT_FLOAT_EQ(test_results[i], fun_results[i]);
}
