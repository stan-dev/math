#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

using namespace Eigen;

TEST(MathFunctions, append_array) {
  std::vector<double> x(3), y(2), z, result;

  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  y[0] = 0.5;
  y[1] = 4.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0]);
  EXPECT_FLOAT_EQ(2.0, result[1]);
  EXPECT_FLOAT_EQ(3.0, result[2]);
  EXPECT_FLOAT_EQ(0.5, result[3]);
  EXPECT_FLOAT_EQ(4.0, result[4]);

  EXPECT_NO_THROW(result = stan::math::append_array(x, z));
  EXPECT_EQ(3, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0]);
  EXPECT_FLOAT_EQ(2.0, result[1]);
  EXPECT_FLOAT_EQ(3.0, result[2]);

  EXPECT_NO_THROW(result = stan::math::append_array(z, z));
  EXPECT_EQ(0, result.size());
}

TEST(MathFunctions, append_array_int_int) {
  std::vector<int> x(2), y(3), result;

  x[0] = 1;
  x[1] = 2;
  y[0] = 3;
  y[1] = 4;
  y[2] = 5;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1, result[0]);
  EXPECT_FLOAT_EQ(2, result[1]);
  EXPECT_FLOAT_EQ(3, result[2]);
  EXPECT_FLOAT_EQ(4, result[3]);
  EXPECT_FLOAT_EQ(5, result[4]);
}

TEST(MathFunctions, append_array_matrix_matrix) {
  std::vector<Matrix<double, Dynamic, Dynamic> > x, y, result, z;

  x.push_back(Matrix<double, Dynamic, Dynamic>(2, 3));
  y.push_back(Matrix<double, Dynamic, Dynamic>(2, 3));
  x[0] << 1.0, 2.0, 3.0,
    4.0, 5.0, 6.0;
  y[0] << 2.0, 3.0, 4.0,
    5.0, 6.0, 7.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(2, result.size());
  EXPECT_FLOAT_EQ(2.0, result[0](0, 1));
  EXPECT_FLOAT_EQ(5.0, result[1](1, 0));

  EXPECT_NO_THROW(result = stan::math::append_array(x, z));
  EXPECT_EQ(1, result.size());

  z.push_back(Matrix<double, Dynamic, Dynamic>(1, 3));
  EXPECT_THROW(result = stan::math::append_array(x, z), std::invalid_argument);

  z.pop_back();
  z.push_back(Matrix<double, Dynamic, Dynamic>(3, 1));
  EXPECT_THROW(result = stan::math::append_array(z, x), std::invalid_argument);
}

TEST(MathFunctions, append_array_vector_vector) {
  std::vector<Matrix<double, Dynamic, 1> > x, y, result, z;

  x.push_back(Matrix<double, Dynamic, 1>(2));
  y.push_back(Matrix<double, Dynamic, 1>(2));
  x[0] << 1.0, 4.0;
  y[0] << 2.0, 5.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(2, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0](0));
  EXPECT_FLOAT_EQ(5.0, result[1](1));

  EXPECT_NO_THROW(result = stan::math::append_array(x, z));
  EXPECT_EQ(1, result.size());

  z.push_back(Matrix<double, Dynamic, 1>(3));
  EXPECT_THROW(result = stan::math::append_array(x, z), std::invalid_argument);
}

TEST(MathFunctions, append_array_row_vector_row_vector) {
  std::vector<Matrix<double, 1, Dynamic> > x, y, result, z;

  x.push_back(Matrix<double, 1, Dynamic>(2));
  y.push_back(Matrix<double, 1, Dynamic>(2));
  x[0] << 1.0, 4.0;
  y[0] << 2.0, 5.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(2, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0](0));
  EXPECT_FLOAT_EQ(5.0, result[1](1));

  EXPECT_NO_THROW(result = stan::math::append_array(x, z));
  EXPECT_EQ(1, result.size());

  z.push_back(Matrix<double, 1, Dynamic>(3));
  EXPECT_THROW(result = stan::math::append_array(x, z), std::invalid_argument);
}

TEST(MathFunctions, append_array_double_int) {
  std::vector<double> x(3), result;
  std::vector<int> y(2);

  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = 3.0;
  y[0] = 5;
  y[1] = 4;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_EQ(5, result.size());
  EXPECT_FLOAT_EQ(1.0, result[0]);
  EXPECT_FLOAT_EQ(2.0, result[1]);
  EXPECT_FLOAT_EQ(3.0, result[2]);
  EXPECT_FLOAT_EQ(5.0, result[3]);
  EXPECT_FLOAT_EQ(4.0, result[4]);
}

TEST(MathFunctions, append_array_nan) {
  double nan = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> x(3), y(2), result;

  x[0] = 1.0;
  x[1] = 2.0;
  x[2] = nan;
  y[0] = 0.5;
  y[1] = 1.0;

  EXPECT_NO_THROW(result = stan::math::append_array(x, y));
  EXPECT_PRED1(boost::math::isnan<double>, result[2]);

  EXPECT_NO_THROW(result = stan::math::append_array(y, x));
  EXPECT_PRED1(boost::math::isnan<double>, result[4]);
}

TEST(MathFunctions, append_array_check_size_vector1) {
  std::vector<std::vector<double> > x(3), y(3), result;

  for(size_t i = 0; i < x.size(); i++) {
    x[i].resize(3);
    y[i].resize(3);
  }
  EXPECT_NO_THROW(result = stan::math::append_array(x, y));

  for(size_t i = 0; i < x.size(); i++) {
    x[i].resize(3);
    y[i].resize(4);
  }
  EXPECT_THROW(result = stan::math::append_array(x, y), std::invalid_argument);
}

TEST(MathFunctions, append_array_check_size_vector2) {
  std::vector<std::vector<std::vector<double> > > x(3), y(3), result;

  for(size_t i = 0; i < 3; i++) {
    x[i].resize(3);
    y[i].resize(3);
    for(size_t j = 0; j < 3; j++) {
      x[i][j].resize(3);
      y[i][j].resize(3);
    }
  }
  EXPECT_NO_THROW(result = stan::math::append_array(x, y));

  for(size_t i = 0; i < 3; i++) {
    x[i].resize(3);
    y[i].resize(3);
    for(size_t j = 0; j < 3; j++) {
      x[i][j].resize(3);
      y[i][j].resize(4);
    }
  }
  EXPECT_THROW(result = stan::math::append_array(x, y), std::invalid_argument);

  for(size_t i = 0; i < 3; i++) {
    x[i].resize(4);
    y[i].resize(3);
    for(size_t j = 0; j < 3; j++) {
      x[i][j].resize(3);
      y[i][j].resize(3);
    }
  }
  EXPECT_THROW(result = stan::math::append_array(x, y), std::invalid_argument);
}
