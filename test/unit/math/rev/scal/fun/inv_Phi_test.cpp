#include <stan/math/rev/scal/fun/inv_Phi.hpp>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>

TEST(AgradRev, inv_Phi) {
  using stan::math::var;

  std::vector<double> p_values;
  p_values.push_back(0.1);
  p_values.push_back(0.5);
  p_values.push_back(0.75);

  // d/dp = sqrt(2 * pi) / exp(normal_log(value_of(x), 0.0, 1.0))
  std::vector<double> dp_values;
  dp_values.push_back(5.69805985611701);
  dp_values.push_back(2.506628274631);
  dp_values.push_back(3.14686508056109);

  for (size_t i = 0; i < p_values.size(); i++) {
    var p, inv_Phi_p;
    AVEC x;
    VEC dp;
    p = p_values[i];
    inv_Phi_p = stan::math::inv_Phi(p);
    x = createAVEC(p);
    inv_Phi_p.grad(x,dp);
    EXPECT_FLOAT_EQ(stan::math::inv_Phi(p.val()), inv_Phi_p.val());
    EXPECT_FLOAT_EQ(dp_values[i], dp[0])
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

