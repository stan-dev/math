#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, dot_product_Eigen) {
  Eigen::Matrix<double, 1, 1> v1;
  v1 << 2.0;
  EXPECT_EQ(4.0, stan::math::dot_product(v1, v1));

  Eigen::Matrix<double, 2, 1> v2;
  v2 << 2.0, 3.0;
  EXPECT_EQ(13.0, stan::math::dot_product(v2, v2));

  Eigen::Matrix<double, 3, 1> v3;
  v3 << 2.0, 3.0, 4.0;
  EXPECT_EQ(29.0, stan::math::dot_product(v3, v3));
}

TEST(MathMatrix, dot_product_Eigen_check) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> v1(1);
  v1 << 1.0;
  Eigen::Matrix<double, Eigen::Dynamic, 1> v2(2);
  v2 << 1.0, 2.0;
  EXPECT_THROW(stan::math::dot_product(v1, v2), std::invalid_argument);
}

TEST(MathMatrix, dot_product_array) {
  std::vector<double> v1 = {2.0};
  EXPECT_EQ(4.0, stan::math::dot_product(v1, v1));

  std::vector<double> v2 = {2.0, 3.0};
  EXPECT_EQ(13.0, stan::math::dot_product(v2, v2));

  std::vector<double> v3 = {2.0, 3.0, 4.0};
  EXPECT_EQ(29.0, stan::math::dot_product(v3, v3));
}

TEST(MathMatrix, dot_product_array_check) {
  stan::math::vector_d v1(2);
  v1 << 1, 2;
  stan::math::vector_d v2(3);
  v2 << 10, 100, 1000;
  EXPECT_THROW(stan::math::dot_product(v1, v2), std::invalid_argument);
}
