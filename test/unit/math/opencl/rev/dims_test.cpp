#ifdef STAN_OPENCL
#include <stan/math/opencl/rev.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MathMatrixRevCL, dimsZero) {
  using stan::math::dims;
  using stan::math::matrix_cl;
  using stan::math::var_value;
  var_value<matrix_cl<double>> m00 = matrix_cl<double>(0, 0);
  std::vector<int> dims00 = dims(m00);
  EXPECT_EQ(2, dims00.size());
  EXPECT_EQ(0, dims00[0]);
  EXPECT_EQ(0, dims00[1]);

  var_value<matrix_cl<double>> m01 = matrix_cl<double>(0, 5);
  std::vector<int> dims01 = dims(m01);
  EXPECT_EQ(2, dims01.size());
  EXPECT_EQ(0, dims01[0]);
  EXPECT_EQ(5, dims01[1]);

  var_value<matrix_cl<double>> m10 = matrix_cl<double>(5, 0);
  std::vector<int> dims10 = dims(m10);
  EXPECT_EQ(2, dims10.size());
  EXPECT_EQ(5, dims10[0]);
  EXPECT_EQ(0, dims10[1]);
}

TEST(MathMatrixRevCL, dimsNonZero) {
  using stan::math::dims;
  using stan::math::matrix_cl;
  using stan::math::var_value;
  var_value<matrix_cl<double>> m1 = matrix_cl<double>(5, 4);
  std::vector<int> dims1 = dims(m1);
  EXPECT_EQ(2, dims1.size());
  EXPECT_EQ(5, dims1[0]);
  EXPECT_EQ(4, dims1[1]);

  var_value<matrix_cl<double>> m2 = matrix_cl<double>(5, 1);
  std::vector<int> dims2 = dims(m2);
  EXPECT_EQ(2, dims2.size());
  EXPECT_EQ(5, dims2[0]);
  EXPECT_EQ(1, dims2[1]);

  var_value<matrix_cl<double>> m3 = matrix_cl<double>(1, 5);
  std::vector<int> dims3 = dims(m3);
  EXPECT_EQ(2, dims3.size());
  EXPECT_EQ(1, dims3[0]);
  EXPECT_EQ(5, dims3[1]);
}
#endif
