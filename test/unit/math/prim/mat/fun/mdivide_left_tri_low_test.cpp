#include <stan/math/prim/mat/fun/mdivide_left_tri_low.hpp>
#include <stan/math/prim/mat/fun/inverse.hpp>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathMatrix, MDivide_LeftTriLow_Values) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a(3, 3);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b(3, 3);
  Eigen::Matrix<double, Eigen::Dynamic, 1> b_vec(3, 1);

  a << 13.5, 0, 0,
       12.5, 3.65, 0,
       1.52, 15.35, 5.69;

  b << 5.56, 2.68, 11.96,
       11.63, 7.52, 1.96,
       6.35, 5.24, 8.65;

  b_vec << 5.56, 11.63, 6.35;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> div_manual(3, 3);
  div_manual << stan::math::inverse(a) * b;

  Eigen::Matrix<double, Eigen::Dynamic, 1> div_vec_manual(3, 1);
  div_vec_manual << stan::math::inverse(a) * b_vec;

  expect_matrix_eq(stan::math::mdivide_left_tri_low(a, b), div_manual);
  expect_matrix_eq(stan::math::mdivide_left_tri_low(a, b_vec), div_vec_manual);
}

TEST(MathMatrix, MDivide_LeftTriLow_Throws) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b_rows_large(6, 3);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b_rows_zero(0, 3);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a(3, 3);
  a << 13.5, 0, 0,
       12.5, 3.65, 0,
       1.52, 15.35, 5.69;

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

  EXPECT_THROW(stan::math::mdivide_left_tri_low(a, b_rows_large),
               std::invalid_argument);
  EXPECT_THROW(stan::math::mdivide_left_tri_low(a, b_rows_zero),
               std::invalid_argument);
  EXPECT_NO_THROW(stan::math::mdivide_left_tri_low(a, b_nan));
  EXPECT_NO_THROW(stan::math::mdivide_left_tri_low(a, b_infty));
  EXPECT_NO_THROW(stan::math::mdivide_left_tri_low(a, b_neg_infty));
}
