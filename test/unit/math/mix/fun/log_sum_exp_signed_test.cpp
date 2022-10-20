#include <test/unit/math/test_ad.hpp>
#include <limits>
#include <vector>

TEST(mathMixScalFun, logSumExp_signed) {
  auto f = [](const int x1_sign, const int x2_sign) {
    return [=](const auto& x1, const auto& x2) {
      stan::return_type_t<decltype(x1), decltype(x2)> ret_val;
      int ret_val_sign;
      std::forward_as_tuple(ret_val, ret_val_sign)
          = stan::math::log_sum_exp_signed(x1, x1_sign, x2, x2_sign);
      return ret_val_sign * stan::math::exp(ret_val);
    };
  };
  std::vector<double> a{0.15, 0.35, 0.51, 0.65, 0.89, 1.0};
  std::vector<double> b{1.4, 1.2, 2.0, 3.0, 3.21, 3.4};

  for (auto&& a_val : a) {
    for (auto&& b_val : b) {
      stan::test::expect_ad(f(1, 1), a_val, b_val);
      stan::test::expect_ad(f(1, -1), a_val, b_val);
      stan::test::expect_ad(f(-1, 1), a_val, b_val);
      stan::test::expect_ad(f(-1, -1), a_val, b_val);
    }
  }
}
