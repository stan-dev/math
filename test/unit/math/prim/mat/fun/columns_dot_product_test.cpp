#include <stan/math/prim/mat/fun/columns_dot_product.hpp>
#include <stan/math/prim/mat/fun/dot_product.hpp>
#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathMatrix, ColumnsDotProduct_Values) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a(3, 3);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b(3, 3);
  Eigen::Matrix<double, Eigen::Dynamic, 1> a1(3, 1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> a2(3, 1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> a3(3, 1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> b1(3, 1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> b2(3, 1);
  Eigen::Matrix<double, Eigen::Dynamic, 1> b3(3, 1);

  a << 13.5, 3.67, 8.98,
       12.5, 3.65, 1.69,
       1.52, 15.35, 5.69;

  a1 << 13.5, 12.5, 1.52;
  a2 << 3.67, 3.65, 15.35;
  a3 << 8.98, 1.69, 5.69;

  b << 5.56, 2.68, 11.96,
       11.63, 7.52, 1.96,
       6.35, 5.24, 8.65;

  b1 << 5.56, 11.63, 6.35;
  b2 << 2.68, 7.52, 5.24;
  b3 << 11.96, 1.96, 8.65;

  Eigen::Matrix<double, 1, Eigen::Dynamic> dot_manual(1, 3);
  dot_manual << (13.5 * 5.56) + (12.5 * 11.63) + (1.52 * 6.35),
                (3.67 * 2.68) + (3.65 * 7.52) + (15.35 * 5.24),
                (8.98 * 11.96) + (1.69 * 1.96) + (5.69 * 8.65);

  Eigen::Matrix<double, 1, Eigen::Dynamic> dot_stan(1, 3);
  dot_stan << stan::math::dot_product(a1, b1),
              stan::math::dot_product(a2, b2),
              stan::math::dot_product(a3, b3);

  Eigen::Matrix<double, 1, Eigen::Dynamic> dot_eigen(1, 3);
  dot_eigen << a.col(0).dot(b.col(0)),
               a.col(1).dot(b.col(1)),
               a.col(2).dot(b.col(2));

  EXPECT_EQ(stan::math::columns_dot_product(a, b), dot_manual);
  EXPECT_EQ(stan::math::columns_dot_product(a, b), dot_stan);
  EXPECT_EQ(stan::math::columns_dot_product(a, b), dot_eigen);
}

TEST(MathMatrix, ColumnsDotProduct_Throws) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a(3, 3);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b_cols_large(3, 6);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b_cols_zero(3, 0);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b_rows_large(6, 3);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b_rows_zero(0, 3);

  double nan = std::numeric_limits<double>::quiet_NaN();
  double infty = std::numeric_limits<double>::infinity();
  double neg_infty = -std::numeric_limits<double>::infinity();

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b_nan(3, 3);
  b_nan << nan, 2.68, 11.96,
           11.63, 7.52, 1.96,
           6.35, 5.24, 8.65;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b_infty(3, 3);
  b_infty << infty, 2.68, 11.96,
             11.63, 7.52, 1.96,
             6.35, 5.24, 8.65;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b_neg_infty(3, 3);
  b_neg_infty << neg_infty, 2.68, 11.96,
                 11.63, 7.52, 1.96,
                 6.35, 5.24, 8.65;

  EXPECT_THROW(stan::math::columns_dot_product(a, b_cols_large),
               std::invalid_argument);
  EXPECT_THROW(stan::math::columns_dot_product(a, b_cols_zero),
               std::invalid_argument);
  EXPECT_THROW(stan::math::columns_dot_product(a, b_rows_large),
               std::invalid_argument);
  EXPECT_THROW(stan::math::columns_dot_product(a, b_rows_zero),
               std::invalid_argument);
  EXPECT_NO_THROW(stan::math::columns_dot_product(a, b_nan));
  EXPECT_NO_THROW(stan::math::columns_dot_product(a, b_infty));
  EXPECT_NO_THROW(stan::math::columns_dot_product(a, b_neg_infty));
}
