#include <stan/math/rev/scal.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/scal/fun/nan_util.hpp>
#include <test/unit/math/rev/scal/util.hpp>
#include <vector>
#include <limits>

TEST(MathFunctions, inv_Phi) {
  using stan::math::Phi;
  using stan::math::inv_Phi;
  using stan::math::var;
  EXPECT_FLOAT_EQ(0.0, inv_Phi(0.5));
  var p = 0.123456789;
  EXPECT_FLOAT_EQ(p.val(), Phi(inv_Phi(p)).val());
  p = 8e-311;
  EXPECT_FLOAT_EQ(p.val(), Phi(inv_Phi(p)).val());
  p = 0.99;
  EXPECT_FLOAT_EQ(p.val(), Phi(inv_Phi(p)).val());

  // breakpoints
  p = 0.02425;
  EXPECT_FLOAT_EQ(p.val(), Phi(inv_Phi(p)).val());
  p = 0.97575;
  EXPECT_FLOAT_EQ(p.val(), Phi(inv_Phi(p)).val());
}
TEST(MathFunctions, inv_Phi_inf) {
  using stan::math::inv_Phi;
  using stan::math::var;
  var p = 7e-311;
  const var inf = std::numeric_limits<var>::infinity();
  EXPECT_EQ(inv_Phi(p), -inf);
  p = 1.0;
  EXPECT_EQ(inv_Phi(p), inf);
}
TEST(MathFunctions, inv_Phi_nan) {
  using stan::math::inv_Phi;
  using stan::math::var;
  var nan = std::numeric_limits<var>::quiet_NaN();
  EXPECT_THROW(inv_Phi(nan), std::domain_error);
  EXPECT_THROW(inv_Phi(-2.0), std::domain_error);
  EXPECT_THROW(inv_Phi(2.0), std::domain_error);
}

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
    y.grad(x, dp);
    EXPECT_FLOAT_EQ(p_values[i], y.val());
    EXPECT_FLOAT_EQ(1.0, dp[0]) << "p = " << p;
  }
}

struct inv_Phi_fun {
  template <typename T0>
  inline T0 operator()(const T0& arg1) const {
    return inv_Phi(arg1);
  }
};

TEST(AgradRev, inv_Phi_NaN) {
  inv_Phi_fun foo;
  test_nan(foo, true, false);
}

TEST(AgradRev, check_varis_on_stack) {
  stan::math::var p = 0.5;
  test::check_varis_on_stack(stan::math::inv_Phi(p));
}
