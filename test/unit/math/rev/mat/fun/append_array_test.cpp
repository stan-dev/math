#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <gtest/gtest.h>

using namespace Eigen;

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

  stan::math::set_zero_all_adjoints();
  result[3].grad();
  EXPECT_FLOAT_EQ(1.0, y[0].adj());

  stan::math::set_zero_all_adjoints();
  result[4].grad();
  EXPECT_FLOAT_EQ(0.0, y[0].adj());
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

  stan::math::set_zero_all_adjoints();
  result[3].grad();
  EXPECT_FLOAT_EQ(0.0, y[0].adj());

  stan::math::set_zero_all_adjoints();
  result[4].grad();
  EXPECT_FLOAT_EQ(0.0, y[0].adj());
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

  stan::math::set_zero_all_adjoints();
  result[4].grad();
  EXPECT_FLOAT_EQ(0.0, y[0].adj());

  stan::math::set_zero_all_adjoints();
  result[1].grad();
  EXPECT_FLOAT_EQ(1.0, x[1].adj());
}

TEST(AgradRev, append_array_matrix_double_matrix_var) {
  std::vector<Matrix<double, Dynamic, Dynamic> > x;
  std::vector<Matrix<stan::math::var, Dynamic, Dynamic> > y, result;

  for(int i = 0; i < 3; i++)
    x.push_back(Matrix<double, Dynamic, Dynamic>::Zero(3, 3));
  for(int i = 0; i < 2; i++)
    y.push_back(Matrix<stan::math::var, Dynamic, Dynamic>::Zero(3, 3));

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

  for(int i = 0; i < 2; i++)
    y[i] = Matrix<stan::math::var, Dynamic, Dynamic>::Zero(2, 2);

  EXPECT_THROW(result = stan::math::append_array(x, y), std::invalid_argument);
}

TEST(AgradRev, append_array_vector_var_vector_var) {
  std::vector<Matrix<stan::math::var, Dynamic, 1> > x, y, result;

  for(int i = 0; i < 3; i++)
    x.push_back(Matrix<stan::math::var, Dynamic, 1>::Zero(3));
  for(int i = 0; i < 2; i++)
    y.push_back(Matrix<stan::math::var, Dynamic, 1>::Zero(3));

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

  for(int i = 0; i < 2; i++)
    y[i] = Matrix<stan::math::var, Dynamic, 1>::Zero(2);

  EXPECT_THROW(result = stan::math::append_array(x, y), std::invalid_argument);
}

TEST(AgradRev, append_array_matrix_types) {
  std::vector<Matrix<double, Dynamic, Dynamic> > xddd(3);
  std::vector<Matrix<stan::math::var, Dynamic, Dynamic> > xvdd(4), rvdd;

  EXPECT_NO_THROW(rvdd = stan::math::append_array(xddd, xvdd));
  EXPECT_NO_THROW(rvdd = stan::math::append_array(xvdd, xddd));
  EXPECT_NO_THROW(rvdd = stan::math::append_array(xvdd, xvdd));

  std::vector<Matrix<double, Dynamic, 1> > xdd1(3);
  std::vector<Matrix<stan::math::var, Dynamic, 1> > xvd1(4), rvd1;

  EXPECT_NO_THROW(rvd1 = stan::math::append_array(xdd1, xvd1));
  EXPECT_NO_THROW(rvd1 = stan::math::append_array(xvd1, xdd1));
  EXPECT_NO_THROW(rvd1 = stan::math::append_array(xvd1, xvd1));

  std::vector<Matrix<double, 1, Dynamic> > xd1d(3);
  std::vector<Matrix<stan::math::var, 1, Dynamic> > xv1d(4), rv1d;

  EXPECT_NO_THROW(rv1d = stan::math::append_array(xd1d, xv1d));
  EXPECT_NO_THROW(rv1d = stan::math::append_array(xv1d, xd1d));
  EXPECT_NO_THROW(rv1d = stan::math::append_array(xv1d, xv1d));
}
