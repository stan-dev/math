#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include <vector>

TEST(MathMatrixPrimMat, dot_product) {
  using stan::math::dot_self;

  Eigen::VectorXd v1(3);
  Eigen::RowVectorXd rv1(3);
  std::vector<double> sv1 = {1.0, 2.0, 3.0};
  v1 << 1.0, 2.0, 3.0;
  rv1 = v1.transpose();

  Eigen::VectorXd v2(3);
  Eigen::RowVectorXd rv2(3);
  std::vector<double> sv2 = {-2.0, -1.0, 3.0};
  v2 << -2.0, -1.0, 3.0;
  rv2 = v2.transpose();

  EXPECT_FLOAT_EQ(5.0, stan::math::dot_product(v1, v2));
  EXPECT_FLOAT_EQ(5.0, stan::math::dot_product(v1, rv2));
  EXPECT_FLOAT_EQ(5.0, stan::math::dot_product(rv1, v2));
  EXPECT_FLOAT_EQ(5.0, stan::math::dot_product(rv1, rv2));
  EXPECT_FLOAT_EQ(5.0, stan::math::dot_product(sv1, sv2));
}

TEST(MathMatrixPrimMat, dot_product_error) {
  using stan::math::dot_self;

  Eigen::VectorXd v1(3);
  Eigen::RowVectorXd rv1(3);
  std::vector<double> sv1 = {1.0, 2.0, 3.0};
  v1 << 1.0, 2.0, 3.0;
  rv1 = v1.transpose();

  Eigen::VectorXd v2(2);
  Eigen::RowVectorXd rv2(2);
  std::vector<double> sv2 = {-2.0, -1.0};
  v2 << -2.0, -1.0;
  rv2 = v2.transpose();

  EXPECT_THROW(stan::math::dot_product(v1, v2), std::invalid_argument);
  EXPECT_THROW(stan::math::dot_product(v1, rv2), std::invalid_argument);
  EXPECT_THROW(stan::math::dot_product(rv1, v2), std::invalid_argument);
  EXPECT_THROW(stan::math::dot_product(rv1, rv2), std::invalid_argument);
  EXPECT_THROW(stan::math::dot_product(sv1, sv2), std::invalid_argument);
}
