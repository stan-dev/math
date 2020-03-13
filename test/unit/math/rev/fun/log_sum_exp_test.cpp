#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/util.hpp>

TEST(log_sum_exp_tests, large_values) {
  using stan::math::var;

  // check combinations of vars and doubles with large argument values
  var a = 1e50;
  var output = stan::math::log_sum_exp(a, a);
  output.grad();
  EXPECT_FLOAT_EQ(output.val(), log(2.0) + value_of(a));
  EXPECT_FLOAT_EQ(a.adj(), 1.0);

  var a2 = 1;
  var a3 = 1e50;
  var output2 = stan::math::log_sum_exp(a2, a3);
  output2.grad();
  EXPECT_FLOAT_EQ(a2.adj(), 0.0);
  EXPECT_FLOAT_EQ(a3.adj(), 1.0);

  var a4 = 1e50;
  var a5 = 1;
  var output3 = stan::math::log_sum_exp(a4, a5);
  output3.grad();
  EXPECT_FLOAT_EQ(a4.adj(), 1.0);
  EXPECT_FLOAT_EQ(a5.adj(), 0.0);


  // check combinations of vars and doubles with large argument values
  var b = 1e20;
  var output4 = stan::math::log_sum_exp(b, b);
  output4.grad();
  EXPECT_FLOAT_EQ(output4.val(), log(2.0) + value_of(b));
  EXPECT_FLOAT_EQ(b.adj(), 1.0);

  var b2 = -2;
  var b3 = 1e20;
  var output5 = stan::math::log_sum_exp(b2, b3);
  output5.grad();
  EXPECT_FLOAT_EQ(b2.adj(), 0.0);
  EXPECT_FLOAT_EQ(b3.adj(), 1.0);

  var b4 = 1e20;
  var b5 = -2;
  var output6 = stan::math::log_sum_exp(b4, b5);
  output6.grad();
  EXPECT_FLOAT_EQ(b4.adj(), 1.0);
  EXPECT_FLOAT_EQ(b5.adj(), 0.0);
}
