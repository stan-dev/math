#include <stan/math/prim/mat/fun/dot_product.hpp>
#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathMatrix, DotProduct_Values) {
  Eigen::Matrix<double, 1, Eigen::Dynamic> r1(1, 3);
  Eigen::Matrix<double, 1, Eigen::Dynamic> r2(1, 3);
  Eigen::Matrix<double, Eigen::Dynamic, 1> c1(3, 1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> c2(3, 1);


  r1 << 13.5, 3.67, 8.98;
  r2 << 12.5, 3.65, 1.69;

  c1 << 5.56, 2.68, 11.96;
  c2 << 11.63, 7.52, 1.96;

  double row_dot_manual = (13.5 * 12.5) + (3.67 * 3.65) + (8.98 * 1.69);
  double col_dot_manual = (5.56 * 11.63) + (2.68 * 7.52) + (11.96 * 1.96);

  double row_dot_eigen = r1.dot(r2);
  double col_dot_eigen = c1.dot(c2);

  EXPECT_EQ(stan::math::dot_product(r1, r2), row_dot_manual);
  EXPECT_EQ(stan::math::dot_product(r1, r2), row_dot_eigen);
  EXPECT_EQ(stan::math::dot_product(c1, c2), col_dot_manual);
  EXPECT_EQ(stan::math::dot_product(c1, c2), col_dot_eigen);
}

TEST(MathMatrix, DotProduct_Throws) {
  Eigen::Matrix<double, Eigen::Dynamic, 1> col_normal(3, 1);
  Eigen::Matrix<double, 1, Eigen::Dynamic> row_normal(1, 3);
  Eigen::Matrix<double, Eigen::Dynamic, 1> col_large(6, 1);
  Eigen::Matrix<double, 1, Eigen::Dynamic> row_large(1, 6);
  Eigen::Matrix<double, Eigen::Dynamic, 1> col_zero(0, 1);
  Eigen::Matrix<double, 1, Eigen::Dynamic> row_zero(1, 0);

  double nan = std::numeric_limits<double>::quiet_NaN();
  double infty = std::numeric_limits<double>::infinity();
  double neg_infty = -std::numeric_limits<double>::infinity();

  Eigen::Matrix<double, Eigen::Dynamic, 1> col_nan(3, 1);
  Eigen::Matrix<double, 1, Eigen::Dynamic> row_nan(1, 3);
  col_nan << nan, 2.68, 11.96;
  row_nan << nan, 2.68, 11.96;

  Eigen::Matrix<double, Eigen::Dynamic, 1> col_infty(3, 1);
  Eigen::Matrix<double, 1, Eigen::Dynamic> row_infty(1, 3);
  col_infty << infty, 2.68, 11.96;
  row_infty << infty, 2.68, 11.96;

  Eigen::Matrix<double, Eigen::Dynamic, 1> col_neg_infty(3, 1);
  Eigen::Matrix<double, 1, Eigen::Dynamic> row_neg_infty(1, 3);
  col_neg_infty << neg_infty, 2.68, 11.96;
  row_neg_infty << neg_infty, 2.68, 11.96;

  EXPECT_THROW(stan::math::dot_product(col_normal, col_large),
               std::invalid_argument);
  EXPECT_THROW(stan::math::dot_product(row_normal, row_large),
               std::invalid_argument);
  EXPECT_THROW(stan::math::dot_product(col_normal, col_zero),
               std::invalid_argument);
  EXPECT_THROW(stan::math::dot_product(row_normal, row_zero),
               std::invalid_argument);
  EXPECT_NO_THROW(stan::math::dot_product(col_normal, col_nan));
  EXPECT_NO_THROW(stan::math::dot_product(row_normal, row_nan));
  EXPECT_NO_THROW(stan::math::dot_product(col_normal, col_infty));
  EXPECT_NO_THROW(stan::math::dot_product(row_normal, row_infty));
  EXPECT_NO_THROW(stan::math::dot_product(col_normal, col_neg_infty));
  EXPECT_NO_THROW(stan::math::dot_product(row_normal, row_neg_infty));
}
