#include <stan/math/prim/scal/fun/selfstart_biexp.hpp>
#include <gtest/gtest.h>
#include <vector>

std::vector<double> test_selfstart_biexp(const std::vector<double>& x,
                                         double A1, double lrc1, double A2,
                                         double lrc2) {
  std::vector<double> value(x.size());
  double exp_lrc1 = std::exp(lrc1);
  double exp_lrc2 = std::exp(lrc2);

  for (std::size_t i = 0; i < x.size(); ++i) {
    value[i] = A1 * std::exp(-exp_lrc1 * x[i]) +
                 A2 * std::exp(-exp_lrc2 * x[i]);
  }
  return value;
}

TEST(MathFunctions, selfstart_biexp) {
  using stan::math::selfstart_biexp;

  double A1 = 3;
  double lrc1 = 1;
  double A2 = 0.6;
  double lrc2 = -1.3;

  std::vector<double> input;
  input.push_back(0);
  input.push_back(-1);
  input.push_back(3.14);

  std::vector<double> test_results =
    test_selfstart_biexp(input, A1, lrc1, A2, lrc2);

  std::vector<double> fun_results =
    selfstart_biexp(input, A1, lrc1, A2, lrc2);

  EXPECT_EQ(test_results.size(), fun_results.size());

  for (std::size_t i = 0; i < test_results.size(); ++i)
    EXPECT_FLOAT_EQ(test_results[i], fun_results[i]);

}
