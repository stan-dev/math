#include <stan/math/fwd/arr.hpp>
#include <gtest/gtest.h>

TEST(AgradFwd, append_array_fvar) {
  std::vector<double> x(3);
  std::vector<stan::math::fvar<double> > y(2), result;

  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  y[0] = 0.5;
  y[0].d_ = 5.0;
  y[1] = 4.0;
  y[1].d_ = 6.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());

  EXPECT_FLOAT_EQ(1.0, result[0].val());
  EXPECT_FLOAT_EQ(2.0, result[1].val());
  EXPECT_FLOAT_EQ(3.0, result[2].val());
  EXPECT_FLOAT_EQ(0.5, result[3].val());
  EXPECT_FLOAT_EQ(4.0, result[4].val());

  EXPECT_FLOAT_EQ(0.0, result[0].tangent());
  EXPECT_FLOAT_EQ(0.0, result[1].tangent());
  EXPECT_FLOAT_EQ(0.0, result[2].tangent());
  EXPECT_FLOAT_EQ(5.0, result[3].tangent());
  EXPECT_FLOAT_EQ(6.0, result[4].tangent());
}
