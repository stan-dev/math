#include <stan/math/rev/mat.hpp>
#include <gtest/gtest.h>
#include <vector>

TEST(MetaTraits, VectorBuilderHelper_false_true) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::VectorBuilderHelper;
  using stan::length;
  using stan::math::var;

  Matrix<var, Dynamic, 1> a_vector(4);
  Matrix<var, 1, Dynamic> a_row_vector(5);

  VectorBuilderHelper<double, false, true> dvv3(length(a_vector));
  EXPECT_THROW(dvv3[0], std::logic_error);
  EXPECT_THROW(dvv3.data(), std::logic_error);

  VectorBuilderHelper<double, false, true> dvv4(length(a_row_vector));
  EXPECT_THROW(dvv4[0], std::logic_error);
  EXPECT_THROW(dvv3.data(), std::logic_error);
}

TEST(MetaTraits, VectorBuilderHelper_true_true) {
  using Eigen::Dynamic;
  using Eigen::Matrix;
  using stan::VectorBuilderHelper;
  using stan::length;
  using stan::math::var;

  Matrix<var, Dynamic, 1> a_vector(4);
  Matrix<var, 1, Dynamic> a_row_vector(5);

  VectorBuilderHelper<double, true, true> dvv3(length(a_vector));
  dvv3[0] = 0.0;
  dvv3[1] = 1.0;
  dvv3[2] = 2.0;
  EXPECT_FLOAT_EQ(0.0, dvv3[0]);
  EXPECT_FLOAT_EQ(1.0, dvv3[1]);
  EXPECT_FLOAT_EQ(2.0, dvv3[2]);
  std::vector<double> data3;
  EXPECT_NO_THROW(data3 = dvv3.data());
  EXPECT_EQ(length(a_vector), data3.size());

  VectorBuilderHelper<double, true, true> dvv4(length(a_row_vector));
  dvv4[0] = 0.0;
  dvv4[1] = 1.0;
  dvv4[2] = 2.0;
  EXPECT_FLOAT_EQ(0.0, dvv4[0]);
  EXPECT_FLOAT_EQ(1.0, dvv4[1]);
  EXPECT_FLOAT_EQ(2.0, dvv4[2]);
  std::vector<double> data4;
  EXPECT_NO_THROW(data4 = dvv4.data());
  EXPECT_EQ(length(a_row_vector), data4.size());
}
