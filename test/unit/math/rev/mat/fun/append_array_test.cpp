#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <gtest/gtest.h>

using namespace Eigen;
using stan::math::var;

TEST(AgradRev, append_array_double_var) {
  std::vector<double> x(3);
  std::vector<var> y(2), result;

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

  stan::math::set_zero_all_adjoints();
  result[3].grad();
  EXPECT_FLOAT_EQ(1.0, y[0].adj());

  stan::math::set_zero_all_adjoints();
  result[4].grad();
  EXPECT_FLOAT_EQ(0.0, y[0].adj());
}

TEST(AgradRev, append_array_var_double) {
  std::vector<double> x(2);
  std::vector<var> y(3), result;

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

  stan::math::set_zero_all_adjoints();
  result[3].grad();
  EXPECT_FLOAT_EQ(0.0, y[0].adj());

  stan::math::set_zero_all_adjoints();
  result[4].grad();
  EXPECT_FLOAT_EQ(0.0, y[0].adj());
}

TEST(AgradRev, append_array_var_var) {
  std::vector<var> x(3), y(2), result;

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

  stan::math::set_zero_all_adjoints();
  result[4].grad();
  EXPECT_FLOAT_EQ(0.0, y[0].adj());

  stan::math::set_zero_all_adjoints();
  result[1].grad();
  EXPECT_FLOAT_EQ(1.0, x[1].adj());
}

TEST(AgradRev, append_array_matrix_double_matrix_var) {
  std::vector<Matrix<double, Dynamic, Dynamic> > x;
  std::vector<Matrix<var, Dynamic, Dynamic> > y, result;

  for (int i = 0; i < 3; i++)
    x.push_back(Matrix<double, Dynamic, Dynamic>::Zero(3, 3));
  for (int i = 0; i < 2; i++)
    y.push_back(Matrix<var, Dynamic, Dynamic>::Zero(3, 3));

  x[0](0, 0) = 1.0;
  y[1](2, 1) = 2.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0](0, 0).val());
  EXPECT_FLOAT_EQ(2.0, result[4](2, 1).val());
  EXPECT_FLOAT_EQ(0.0, result[4](2, 2).val());

  stan::math::set_zero_all_adjoints();
  result[4](2, 1).grad();
  EXPECT_FLOAT_EQ(0.0, y[1](1, 2).adj());
  EXPECT_FLOAT_EQ(1.0, y[1](2, 1).adj());

  for (int i = 0; i < 2; i++)
    y[i] = Matrix<var, Dynamic, Dynamic>::Zero(2, 2);

  EXPECT_THROW(result = stan::math::append_array(x, y), std::invalid_argument);
}

TEST(AgradRev, append_array_vector_var_vector_var) {
  std::vector<Matrix<var, Dynamic, 1> > x, y, result;

  for (int i = 0; i < 3; i++)
    x.push_back(Matrix<var, Dynamic, 1>::Zero(3));
  for (int i = 0; i < 2; i++)
    y.push_back(Matrix<var, Dynamic, 1>::Zero(3));

  x[0](0) = 1.0;
  y[1](2) = 2.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0](0).val());
  EXPECT_FLOAT_EQ(0.0, result[4](1).val());
  EXPECT_FLOAT_EQ(2.0, result[4](2).val());

  stan::math::set_zero_all_adjoints();
  result[4](2).grad();
  EXPECT_FLOAT_EQ(0.0, y[1](1).adj());
  EXPECT_FLOAT_EQ(1.0, y[1](2).adj());

  for (int i = 0; i < 2; i++)
    y[i] = Matrix<var, Dynamic, 1>::Zero(2);

  EXPECT_THROW(result = stan::math::append_array(x, y), std::invalid_argument);
}
