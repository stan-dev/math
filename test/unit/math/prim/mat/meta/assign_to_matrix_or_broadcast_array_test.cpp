#include <gtest/gtest.h>
#include <stan/math/prim/mat.hpp>
#include <stan/math/prim/scal.hpp>

TEST(MetaTraits, AssignToMatrixOrBroadcastArrayTestBA) {
  double two = 2;
  stan::math::internal::broadcast_array<double> BA(two);
  Eigen::Matrix<double, -1, 1> EM(3, 1);
  EM << 1, 2, 3;
  stan::math::assign_to_matrix_or_broadcast_array(BA, EM);
  EXPECT_FLOAT_EQ(BA[0], 1);
  EXPECT_FLOAT_EQ(EM[0], 1);
  EXPECT_FLOAT_EQ(EM[1], 2);
  EXPECT_FLOAT_EQ(EM[2], 3);
}

TEST(MetaTraits, AssignToMatrixOrBroadcastArrayTestEM) {
  Eigen::Matrix<double, -1, 1> EM(3, 1);
  EM << 1, 2, 3;
  Eigen::Matrix<double, -1, 1> EM2(3, 1);
  EM2 << 4, 5, 6;
  stan::math::assign_to_matrix_or_broadcast_array(EM, EM2);
  EXPECT_FLOAT_EQ(EM[0], 4);
  EXPECT_FLOAT_EQ(EM[1], 5);
  EXPECT_FLOAT_EQ(EM[2], 6);
  EXPECT_FLOAT_EQ(EM2[0], 4);
  EXPECT_FLOAT_EQ(EM2[1], 5);
  EXPECT_FLOAT_EQ(EM2[2], 6);
}
