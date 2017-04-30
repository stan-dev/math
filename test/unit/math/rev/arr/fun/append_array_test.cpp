#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, append_array_var) {
  std::vector<double> x(3), dr;
  std::vector<stan::math::var> y(2), result;

  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  y[0] = 0.5;
  y[1] = 4.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0].val());
  EXPECT_FLOAT_EQ(2.0, result[1].val());
  EXPECT_FLOAT_EQ(3.0, result[2].val());
  EXPECT_FLOAT_EQ(0.5, result[3].val());
  EXPECT_FLOAT_EQ(4.0, result[4].val());

  result[0].grad(y, dr);
  EXPECT_FLOAT_EQ(0.0, dr[0]);
  EXPECT_FLOAT_EQ(0.0, dr[1]);
  stan::math::set_zero_all_adjoints();
  result[1].grad(y, dr);
  EXPECT_FLOAT_EQ(0.0, dr[0]);
  EXPECT_FLOAT_EQ(0.0, dr[1]);
  stan::math::set_zero_all_adjoints();
  result[2].grad(y, dr);
  EXPECT_FLOAT_EQ(0.0, dr[0]);
  EXPECT_FLOAT_EQ(0.0, dr[1]);
  stan::math::set_zero_all_adjoints();
  result[3].grad(y, dr);
  EXPECT_FLOAT_EQ(1.0, dr[0]);
  EXPECT_FLOAT_EQ(0.0, dr[1]);
  stan::math::set_zero_all_adjoints();
  result[4].grad(y, dr);
  EXPECT_FLOAT_EQ(0.0, dr[0]);
  EXPECT_FLOAT_EQ(1.0, dr[1]);
}
