#include <stan/math/prim/mat/fun/rows_dot_product.hpp>
#include <stan/math/prim/mat/fun/dot_product.hpp>
#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathMatrix, RowsDotProduct_Values) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a(3, 3);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b(3, 3);
  Eigen::Matrix<double, 1, Eigen::Dynamic> a1(1, 3);
  Eigen::Matrix<double, 1, Eigen::Dynamic> a2(1, 3);
  Eigen::Matrix<double, 1, Eigen::Dynamic> a3(1, 3);
  Eigen::Matrix<double, 1, Eigen::Dynamic> b1(1, 3);
  Eigen::Matrix<double, 1, Eigen::Dynamic> b2(1, 3);
  Eigen::Matrix<double, 1, Eigen::Dynamic> b3(1, 3);

  a << 13.5, 3.67, 8.98,
       12.5, 3.65, 1.69,
       1.52, 15.35, 5.69;

  a1 << 13.5, 3.67, 8.98;
  a2 << 12.5, 3.65, 1.69;
  a3 << 1.52, 15.35, 5.69;

  b << 5.56, 2.68, 11.96,
       11.63, 7.52, 1.96,
       6.35, 5.24, 8.65;

  b1 << 5.56, 2.68, 11.96;
  b2 << 11.63, 7.52, 1.96;
  b3 << 6.35, 5.24, 8.65;

  Eigen::Matrix<double, Eigen::Dynamic,  1> dot_manual(3, 1);
  dot_manual << (13.5 * 5.56) + (3.67 * 2.68) + (8.98 * 11.96),
                (12.5 * 11.63) + (3.65 * 7.52) + (1.69 * 1.96),
                (1.52 * 6.35) + (15.35 * 5.24) + (5.69 * 8.65);

  Eigen::Matrix<double, Eigen::Dynamic, 1> dot_stan(3, 1);
  dot_stan << stan::math::dot_product(a1, b1),
              stan::math::dot_product(a2, b2),
              stan::math::dot_product(a3, b3);

  Eigen::Matrix<double, Eigen::Dynamic, 1> dot_eigen(3, 1);
  dot_eigen << a.row(0).dot(b.row(0)),
               a.row(1).dot(b.row(1)),
               a.row(2).dot(b.row(2));

  EXPECT_EQ(stan::math::rows_dot_product(a, b), dot_manual);
  EXPECT_EQ(stan::math::rows_dot_product(a, b), dot_stan);
  EXPECT_EQ(stan::math::rows_dot_product(a, b), dot_eigen);
}

TEST(MathMatrix, RowsDotProduct_Throws) {
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

  EXPECT_THROW(stan::math::rows_dot_product(a, b_cols_large),
               std::invalid_argument);
  EXPECT_THROW(stan::math::rows_dot_product(a, b_cols_zero),
               std::invalid_argument);
  EXPECT_THROW(stan::math::rows_dot_product(a, b_rows_large),
               std::invalid_argument);
  EXPECT_THROW(stan::math::rows_dot_product(a, b_rows_zero),
               std::invalid_argument);
  EXPECT_NO_THROW(stan::math::rows_dot_product(a, b_nan));
  EXPECT_NO_THROW(stan::math::rows_dot_product(a, b_infty));
  EXPECT_NO_THROW(stan::math::rows_dot_product(a, b_neg_infty));
}
