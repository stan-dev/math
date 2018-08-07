#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/scal.hpp>
#include <vector>

TEST(MetaTraits, DuplicateIfScalarTestSca) {
  double two = 2;
  int length = 3;
  Eigen::Matrix<double, Eigen::Dynamic, 1> output
    = stan::math::duplicate_if_scalar(two, length);
  EXPECT_FLOAT_EQ(output.size(), length);
  for (size_t ii = 0; ii < length; ii++) {
    EXPECT_FLOAT_EQ(output[ii], two);
  }
}

TEST(MetaTraits, DuplicateIfScalarTestVec) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> vector(3, 1);
  vector << 1, 2, 3;
  Eigen::Matrix<double, Eigen::Dynamic, 1> output
    = stan::math::duplicate_if_scalar(vector, 3);
  EXPECT_FLOAT_EQ(output.size(), 3);
  for (size_t ii = 0; ii < 3; ii++) {
    EXPECT_FLOAT_EQ(output[ii], vector[ii]);
  }
}
