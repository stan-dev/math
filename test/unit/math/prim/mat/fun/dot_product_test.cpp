#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, dot_product_Eigen) {
  Eigen::Matrix<double, 1, 1> v1;
  v1 << 2.0;
  EXPECT_NEAR(4.0, stan::math::dot_product(v1, v1), 1E-12);

  Eigen::Matrix<double, 2, 1> v2;
  v2 << 2.0, 3.0;
  EXPECT_NEAR(13.0, stan::math::dot_product(v2, v2), 1E-12);

  Eigen::Matrix<double, 3, 1> v3;
  v3 << 2.0, 3.0, 4.0;
  EXPECT_NEAR(29.0, stan::math::dot_product(v3, v3), 1E-12);
}

TEST(MathMatrix, dot_product_array) {
  std::vector<double> v1 = {2.0};
  EXPECT_NEAR(4.0, stan::math::dot_product(v1, v1), 1E-12);

  std::vector<double> v2 = {2.0, 3.0};
  EXPECT_NEAR(13.0, stan::math::dot_product(v2, v2), 1E-12);

  std::vector<double> v3 = {2.0, 3.0, 4.0};
  EXPECT_NEAR(29.0, stan::math::dot_product(v3, v3), 1E-12);
}
