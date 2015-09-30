#include <stan/math/prim/scal/fun/selfstart_log.hpp>
#include <gtest/gtest.h>
#include <vector>

std::vector<double> test_selfstart_four_param_log(
                                         const std::vector<double>& x,
                                         double A, double B, double xmid,
                                         double scal) {
  std::vector<double> value(x.size());

  for (std::size_t i = 0; i < x.size(); ++i) {
    value[i] = A + (B - A)/(1 + std::exp((xmid - x[i])/scal));
  }
  return value;
}

std::vector<double> test_selfstart_log(const std::vector<double>& x,
                                       double Asym, double xmid, double scal) {
  std::vector<double> value(x.size());

  for (std::size_t i = 0; i < x.size(); ++i) {
    value[i] = Asym/(1 + std::exp((xmid - x[i])/scal));
  }
  return value;
}


TEST(MathFunctions, selfstart_four_param_log) {
  using stan::math::selfstart_four_param_log;

  double A = 13;
  double B = 368;
  double xmid = 14;
  double scal = 6;

  std::vector<double> input;
  input.push_back(0);
  input.push_back(-1);
  input.push_back(3.14);

  std::vector<double> test_results =
    test_selfstart_four_param_log(input, A, B, xmid, scal);

  std::vector<double> fun_results =
    selfstart_four_param_log(input, A, B, xmid, scal);

  EXPECT_EQ(test_results.size(), fun_results.size());

  for (std::size_t i = 0; i < test_results.size(); ++i)
    EXPECT_FLOAT_EQ(test_results[i], fun_results[i]);
}

TEST(MathFunctions, selfstart_log) {
  using stan::math::selfstart_log;

  double Asym = 368;
  double xmid = 14;
  double scal = 6;

  std::vector<double> input;
  input.push_back(0);
  input.push_back(-1);
  input.push_back(3.14);

  std::vector<double> test_results =
    test_selfstart_log(input, Asym, xmid, scal);

  std::vector<double> fun_results =
    selfstart_log(input, Asym, xmid, scal);

  EXPECT_EQ(test_results.size(), fun_results.size());

  for (std::size_t i = 0; i < test_results.size(); ++i)
    EXPECT_FLOAT_EQ(test_results[i], fun_results[i]);
}

