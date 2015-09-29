#include <stan/math/prim/scal/fun/selfstart_first_order_compart.hpp>
#include <gtest/gtest.h>
#include <vector>

std::vector<double> test_selfstart_first_order_compart(
                                         const std::vector<double>& x,
                                         double Dose, double lKe, double lKa,
                                         double lCl) {
  std::vector<double> value(x.size());
  double exp_lke = std::exp(lKe);
  double exp_lka = std::exp(lKa);
  double dose_const = Dose * exp(lKe + lKa - lCl)/(exp_lka - exp_lke);

  for (std::size_t i = 0; i < x.size(); ++i) {
    value[i] = dose_const *
                 (std::exp(-exp_lke * x[i]) - std::exp(-exp_lka * x[i]));
  }
  return value;
}

TEST(MathFunctions, selfstart_first_order_compart) {
  using stan::math::selfstart_first_order_compart;

  double dose = 500;
  double lKe = -2.5;
  double lKa = 0.5;
  double lCl = -3;

  std::vector<double> input;
  input.push_back(0);
  input.push_back(-1);
  input.push_back(3.14);

  std::vector<double> test_results =
    test_selfstart_first_order_compart(input, dose, lKe, lKa, lCl);

  std::vector<double> fun_results =
    selfstart_first_order_compart(input, dose, lKe, lKa, lCl);

  EXPECT_EQ(test_results.size(), fun_results.size());

  for (std::size_t i = 0; i < test_results.size(); ++i)
    EXPECT_FLOAT_EQ(test_results[i], fun_results[i]);
}
