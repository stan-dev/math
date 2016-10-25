#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>

TEST(AgradRevMatrix, divide_scalar) {
  using stan::math::divide;
  double d1, d2;
  AVAR   v1, v2;

  d1 = 10;
  v1 = 10;
  d2 = -2;
  v2 = -2;
  
  EXPECT_FLOAT_EQ(-5, divide(d1, d2));
  EXPECT_FLOAT_EQ(-5, divide(d1, v2).val());
  EXPECT_FLOAT_EQ(-5, divide(v1, d2).val());
  EXPECT_FLOAT_EQ(-5, divide(v1, v2).val());

  d2 = 0;
  v2 = 0;

  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), divide(d1, d2));
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), divide(d1, v2).val());
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), divide(v1, d2).val());
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), divide(v1, v2).val());

  d1 = 0;
  v1 = 0;
  EXPECT_TRUE(std::isnan(divide(d1, d2)));
  EXPECT_TRUE(std::isnan(divide(d1, v2).val()));
  EXPECT_TRUE(std::isnan(divide(v1, d2).val()));
  EXPECT_TRUE(std::isnan(divide(v1, v2).val()));
}
TEST(AgradRevMatrix, divide_vector) {
  using stan::math::divide;
  using stan::math::vector_d;
  using stan::math::vector_v;

  vector_d d1(3);
  vector_v v1(3);
  double d2;
  AVAR v2;
  
  d1 << 100, 0, -3;
  v1 << 100, 0, -3;
  d2 = -2;
  v2 = -2;
  
  vector_d output_d;
  output_d = divide(d1, d2);
  EXPECT_FLOAT_EQ(-50, output_d(0));
  EXPECT_FLOAT_EQ(  0, output_d(1));
  EXPECT_FLOAT_EQ(1.5, output_d(2));

  vector_v output;
  output = divide(d1, v2);
  EXPECT_FLOAT_EQ(-50, output(0).val());
  EXPECT_FLOAT_EQ(  0, output(1).val());
  EXPECT_FLOAT_EQ(1.5, output(2).val());

  output = divide(v1, d2);
  EXPECT_FLOAT_EQ(-50, output(0).val());
  EXPECT_FLOAT_EQ(  0, output(1).val());
  EXPECT_FLOAT_EQ(1.5, output(2).val());

  output = divide(v1, v2);
  EXPECT_FLOAT_EQ(-50, output(0).val());
  EXPECT_FLOAT_EQ(  0, output(1).val());
  EXPECT_FLOAT_EQ(1.5, output(2).val());


  d2 = 0;
  v2 = 0;
  output_d = divide(d1, d2);
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output_d(0));
  EXPECT_TRUE (std::isnan(output_d(1)));
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), output_d(2));

  output = divide(d1, v2);
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output(0).val());
  EXPECT_TRUE (std::isnan(output(1).val()));
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), output(2).val());

  output = divide(v1, d2);
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output(0).val());
  EXPECT_TRUE (std::isnan(output(1).val()));
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), output(2).val());

  output = divide(v1, v2);
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output(0).val());
  EXPECT_TRUE (std::isnan(output(1).val()));
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), output(2).val());
}
TEST(AgradRevMatrix, divide_rowvector) {
  using stan::math::divide;
  using stan::math::row_vector_d;
  using stan::math::row_vector_v;

  row_vector_d d1(3);
  row_vector_v v1(3);
  double d2;
  AVAR v2;
  
  d1 << 100, 0, -3;
  v1 << 100, 0, -3;
  d2 = -2;
  v2 = -2;
  
  row_vector_d output_d = divide(d1, d2);
  EXPECT_FLOAT_EQ(-50, output_d(0));
  EXPECT_FLOAT_EQ(  0, output_d(1));
  EXPECT_FLOAT_EQ(1.5, output_d(2));

  row_vector_v output;
  output = divide(d1, v2);
  EXPECT_FLOAT_EQ(-50, output(0).val());
  EXPECT_FLOAT_EQ(  0, output(1).val());
  EXPECT_FLOAT_EQ(1.5, output(2).val());

  output = divide(v1, d2);
  EXPECT_FLOAT_EQ(-50, output(0).val());
  EXPECT_FLOAT_EQ(  0, output(1).val());
  EXPECT_FLOAT_EQ(1.5, output(2).val());

  output = divide(v1, v2);
  EXPECT_FLOAT_EQ(-50, output(0).val());
  EXPECT_FLOAT_EQ(  0, output(1).val());
  EXPECT_FLOAT_EQ(1.5, output(2).val());

  d2 = 0;
  v2 = 0;
  output_d = divide(d1, d2);
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output_d(0));
  EXPECT_TRUE(std::isnan(output_d(1)));
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), output_d(2));

  output = divide(d1, v2);
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output(0).val());
  EXPECT_TRUE(std::isnan(output(1).val()));
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), output(2).val());

  output = divide(v1, d2);
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output(0).val());
  EXPECT_TRUE (std::isnan(output(1).val()));
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), output(2).val());

  output = divide(v1, v2);
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output(0).val());
  EXPECT_TRUE (std::isnan(output(1).val()));
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), output(2).val());
}
TEST(AgradRevMatrix, divide_matrix) {
  using stan::math::divide;
  using stan::math::matrix_d;
  using stan::math::matrix_v;

  matrix_d d1(2,2);
  matrix_v v1(2,2);
  double d2;
  AVAR v2;
  
  d1 << 100, 0, -3, 4;
  v1 << 100, 0, -3, 4;
  d2 = -2;
  v2 = -2;
  
  matrix_d output_d = divide(d1, d2);
  EXPECT_FLOAT_EQ(-50, output_d(0,0));
  EXPECT_FLOAT_EQ(  0, output_d(0,1));
  EXPECT_FLOAT_EQ(1.5, output_d(1,0));
  EXPECT_FLOAT_EQ( -2, output_d(1,1));

  matrix_v output;
  output = divide(d1, v2);
  EXPECT_FLOAT_EQ(-50, output(0,0).val());
  EXPECT_FLOAT_EQ(  0, output(0,1).val());
  EXPECT_FLOAT_EQ(1.5, output(1,0).val());
  EXPECT_FLOAT_EQ( -2, output(1,1).val());
  
  output = divide(v1, d2);
  EXPECT_FLOAT_EQ(-50, output(0,0).val());
  EXPECT_FLOAT_EQ(  0, output(0,1).val());
  EXPECT_FLOAT_EQ(1.5, output(1,0).val());
  EXPECT_FLOAT_EQ( -2, output(1,1).val());
  
  output = divide(v1, v2);
  EXPECT_FLOAT_EQ(-50, output(0,0).val());
  EXPECT_FLOAT_EQ(  0, output(0,1).val());
  EXPECT_FLOAT_EQ(1.5, output(1,0).val());
  EXPECT_FLOAT_EQ( -2, output(1,1).val());

  d2 = 0;
  v2 = 0;
  output_d = divide(d1, d2);
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output_d(0,0));
  EXPECT_TRUE(std::isnan(output_d(0,1)));
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), output_d(1,0));
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output_d(1,1));

  output = divide(d1, v2);
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output(0,0).val());
  EXPECT_TRUE (std::isnan(output(0,1).val()));
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), output(1,0).val());
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output(1,1).val());

  output = divide(v1, d2);
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output(0,0).val());
  EXPECT_TRUE (std::isnan(output(0,1).val()));
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), output(1,0).val());
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output(1,1).val());

  output = divide(v1, v2);
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output(0,0).val());
  EXPECT_TRUE (std::isnan(output(0,1).val()));
  EXPECT_FLOAT_EQ(-std::numeric_limits<double>::infinity(), output(1,0).val());
  EXPECT_FLOAT_EQ(std::numeric_limits<double>::infinity(), output(1,1).val());
}

TEST(AgradRevMatrix, check_varis_on_stack) {
  using stan::math::value_of;
  stan::math::var x = 10;
  stan::math::vector_v v(3);
  v << -100, 0, 1;
  stan::math::row_vector_v rv(3);
  rv << -100, 0, 1;
  stan::math::matrix_v m(2, 3);
  m << -100, 0, 1, 20, -40, 2;
  
  test::check_varis_on_stack(stan::math::divide(v, x));
  test::check_varis_on_stack(stan::math::divide(v, value_of(x)));
  test::check_varis_on_stack(stan::math::divide(value_of(v), x));

  test::check_varis_on_stack(stan::math::divide(rv, x));
  test::check_varis_on_stack(stan::math::divide(rv, value_of(x)));
  test::check_varis_on_stack(stan::math::divide(value_of(rv), x));

  test::check_varis_on_stack(stan::math::divide(m, x));
  test::check_varis_on_stack(stan::math::divide(m, value_of(x)));
  test::check_varis_on_stack(stan::math::divide(value_of(m), x));
}
