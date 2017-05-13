#include <stan/math/rev/arr.hpp>
#include <gtest/gtest.h>

TEST(AgradRev, append_array_double_var) {
  std::vector<double> x(3);
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

  for(size_t i = 0; i < result.size(); i++) {
    stan::math::set_zero_all_adjoints();

    result[i].grad();

    if(i == 3)
      EXPECT_FLOAT_EQ(1.0, y[0].adj());
    else
      EXPECT_FLOAT_EQ(0.0, y[0].adj());

    if(i == 4)
      EXPECT_FLOAT_EQ(1.0, y[1].adj());
    else
      EXPECT_FLOAT_EQ(0.0, y[1].adj());
  }
}

TEST(AgradRev, append_array_var_double) {
  std::vector<double> x(2);
  std::vector<stan::math::var> y(3), result;

  x[0] = 1.0;
  x[1] = 2.0;
  y[0] = 5.0;
  y[1] = 6.0;
  y[2] = 7.0;

  EXPECT_NO_THROW(result = stan::math::append_array(y, x));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(5.0, result[0].val());
  EXPECT_FLOAT_EQ(6.0, result[1].val());
  EXPECT_FLOAT_EQ(7.0, result[2].val());
  EXPECT_FLOAT_EQ(1.0, result[3].val());
  EXPECT_FLOAT_EQ(2.0, result[4].val());

  for(size_t i = 0; i < result.size(); i++) {
    stan::math::set_zero_all_adjoints();

    result[i].grad();

    if(i == 0)
      EXPECT_FLOAT_EQ(1.0, y[0].adj());
    else
      EXPECT_FLOAT_EQ(0.0, y[0].adj());

    if(i == 1)
      EXPECT_FLOAT_EQ(1.0, y[1].adj());
    else
      EXPECT_FLOAT_EQ(0.0, y[1].adj());

    if(i == 2)
      EXPECT_FLOAT_EQ(1.0, y[2].adj());
    else
      EXPECT_FLOAT_EQ(0.0, y[2].adj());
  }
}

TEST(AgradRev, append_array_var_var) {
  std::vector<stan::math::var> x(3), y(2), result;

  x[0] = 5.0;
  x[1] = 6.0;
  x[2] = 7.0;
  y[0] = 0.5;
  y[1] = 4.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(5.0, result[0].val());
  EXPECT_FLOAT_EQ(6.0, result[1].val());
  EXPECT_FLOAT_EQ(7.0, result[2].val());
  EXPECT_FLOAT_EQ(0.5, result[3].val());
  EXPECT_FLOAT_EQ(4.0, result[4].val());

  for(size_t i = 0; i < result.size(); i++) {
    stan::math::set_zero_all_adjoints();

    result[i].grad();

    if(i == 0)
      EXPECT_FLOAT_EQ(1.0, x[0].adj());
    else
      EXPECT_FLOAT_EQ(0.0, x[0].adj());

    if(i == 1)
      EXPECT_FLOAT_EQ(1.0, x[1].adj());
    else
      EXPECT_FLOAT_EQ(0.0, x[1].adj());

    if(i == 2)
      EXPECT_FLOAT_EQ(1.0, x[2].adj());
    else
      EXPECT_FLOAT_EQ(0.0, x[2].adj());

    if(i == 3)
      EXPECT_FLOAT_EQ(1.0, y[0].adj());
    else
      EXPECT_FLOAT_EQ(0.0, y[0].adj());

    if(i == 4)
      EXPECT_FLOAT_EQ(1.0, y[1].adj());
    else
      EXPECT_FLOAT_EQ(0.0, y[1].adj());
  }
}
