#include <stan/math/rev/scal/fun/inv_Phi.hpp>
#include <stan/math/rev/scal/fun/Phi.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>

TEST(AgradRev, inv_Phi) {
  using stan::math::var;

  std::vector<double> p_values;
  p_values.push_back(0.1);
  p_values.push_back(0.5);
  p_values.push_back(0.75);

  for (size_t i = 0; i < p_values.size(); i++) {
    var p, y;
    AVEC x;
    VEC dp;
    p = p_values[i];
    y = stan::math::Phi(stan::math::inv_Phi(p));
    x = createAVEC(p);
    y.grad(x,dp);
    EXPECT_FLOAT_EQ(p_values[i], y.val());
    EXPECT_FLOAT_EQ(1.0, dp[0])
      << "p = " << p;
  }
}

struct inv_Phi_fun {
  template <typename T0>
  inline T0
  operator()(const T0& arg1) const {
    return inv_Phi(arg1);
  }
};

TEST(AgradRev,inv_Phi_NaN) {
  inv_Phi_fun inv_Phi_;
  test_nan(inv_Phi_,true,false);
}

