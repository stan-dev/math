#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixPrimMat, mean) {
  using stan::math::mean;
  std::vector<double> x;
  EXPECT_THROW(mean(x), std::invalid_argument);
  x.push_back(1.0);
  EXPECT_FLOAT_EQ(1.0, mean(x));
  x.push_back(2.0);
  EXPECT_NEAR(1.5, mean(x), 0.000001);
  x.push_back(3.0);
  EXPECT_FLOAT_EQ(2.0, mean(x));

  stan::math::vector_d v;
  EXPECT_THROW(mean(v), std::invalid_argument);
  v = stan::math::vector_d(1);
  v << 1.0;
  EXPECT_FLOAT_EQ(1.0, mean(v));
  v = stan::math::vector_d(2);
  v << 1.0, 2.0;
  EXPECT_NEAR(1.5, mean(v), 0.000001);
  v = stan::math::vector_d(3);
  v << 1.0, 2.0, 3.0;
  EXPECT_FLOAT_EQ(2.0, mean(v));

  stan::math::row_vector_d rv;
  EXPECT_THROW(mean(rv), std::invalid_argument);
  rv = stan::math::row_vector_d(1);
  rv << 1.0;
  EXPECT_FLOAT_EQ(1.0, mean(rv));
  rv = stan::math::row_vector_d(2);
  rv << 1.0, 2.0;
  EXPECT_NEAR(1.5, mean(rv), 0.000001);
  rv = stan::math::row_vector_d(3);
  rv << 1.0, 2.0, 3.0;
  EXPECT_FLOAT_EQ(2.0, mean(rv));

  stan::math::matrix_d m;
  EXPECT_THROW(mean(m), std::invalid_argument);
  m = stan::math::matrix_d(1, 1);
  m << 1.0;
  EXPECT_FLOAT_EQ(1.0, mean(m));
  m = stan::math::matrix_d(2, 3);
  m << 1.0, 2.0, 4.0, 9.0, 16.0, 25.0;
  EXPECT_FLOAT_EQ(9.5, mean(m));
}

TEST(MathMatrixPrimMat, mean_exception) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  Matrix<double, Dynamic, Dynamic> m;
  Matrix<double, Dynamic, 1> v;
  Matrix<double, 1, Dynamic> rv;

  EXPECT_THROW(stan::math::mean(m), std::invalid_argument);
  EXPECT_THROW(stan::math::mean(v), std::invalid_argument);
  EXPECT_THROW(stan::math::mean(rv), std::invalid_argument);

  Matrix<double, Dynamic, Dynamic> m_nz(2, 3);
  Matrix<double, Dynamic, 1> v_nz(2);
  Matrix<double, 1, Dynamic> rv_nz(3);

  EXPECT_NO_THROW(stan::math::mean(m_nz));
  EXPECT_NO_THROW(stan::math::mean(v_nz));
  EXPECT_NO_THROW(stan::math::mean(rv_nz));
}
