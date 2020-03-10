#include <stan/math/rev.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/fun/util.hpp>
#include <test/unit/math/rev/util.hpp>

TEST(log_sum_exp_tests, large_values) {
 using stan::math::var;

 var a = 1e50;
 var output = stan::math::log_sum_exp(a, a);
 output.grad();
 EXPECT_FLOAT_EQ(output.val(), log(2.0) + value_of(a));
 EXPECT_FLOAT_EQ(a.adj(), 1.0);

 var a1 = 1e50;
 var a2 = 1;
 var output2 = stan::math::log_sum_exp(a1, a2);
 output2.grad();
 EXPECT_FLOAT_EQ(a1.adj(), 1.0);
 EXPECT_FLOAT_EQ(a2.adj(), 0.0);

 var a3 = 1;
 var a4 = 1e50;
 var output3 = stan::math::log_sum_exp(a3, a4);
 output3.grad();
 EXPECT_FLOAT_EQ(a3.adj(), 0.0);
 EXPECT_FLOAT_EQ(a4.adj(), 1.0);
}
