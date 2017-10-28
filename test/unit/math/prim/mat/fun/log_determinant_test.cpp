#include <stan/math/prim/mat/fun/log_determinant.hpp>
#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <limits>

TEST(MathMatrix, LogDeterminant_Values) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a(3, 3);
  a << 1.52, 1.69, 11.98,
       12.56, 37.56, 31.69,
       15.42, 43.31, 5.69;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> b(3, 3);
  b << 5.56, 2.68, 11.96,
       11.63, 7.52, 1.96,
       6.35, 5.24, 8.65;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a_times_b(3, 3);
  a_times_b << a * b;

  double det_manual = (1.52 * 37.56 * 5.69) + (1.69 * 31.69 * 15.42)
                       + (11.98 * 12.56 * 43.31) - (11.98 * 37.56 * 15.42)
                       - (1.69 * 12.56 * 5.69)  - (1.52 * 31.69 * 43.31);
  double sum_logdet = stan::math::log_determinant(a)
                      + stan::math::log_determinant(b);

  EXPECT_FLOAT_EQ(stan::math::log_determinant(a), log(abs(det_manual)));
  EXPECT_FLOAT_EQ(stan::math::log_determinant(a_times_b), sum_logdet);
}

TEST(MathMatrix, LogDeterminant_Throws) {
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a_cols_large(3, 6);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a_cols_zero(3, 0);

  double nan = std::numeric_limits<double>::quiet_NaN();
  double infty = std::numeric_limits<double>::infinity();
  double neg_infty = -std::numeric_limits<double>::infinity();

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a_nan(3, 3);
  a_nan << nan, 2.68, 11.96,
           11.63, 7.52, 1.96,
           6.35, 5.24, 8.65;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a_infty(3, 3);
  a_infty << infty, 2.68, 11.96,
             11.63, 7.52, 1.96,
             6.35, 5.24, 8.65;

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> a_neg_infty(3, 3);
  a_neg_infty << neg_infty, 2.68, 11.96,
                 11.63, 7.52, 1.96,
                 6.35, 5.24, 8.65;

  EXPECT_THROW(stan::math::log_determinant(a_cols_large),
               std::invalid_argument);
  EXPECT_THROW(stan::math::log_determinant(a_cols_zero),
               std::invalid_argument);
  EXPECT_NO_THROW(stan::math::log_determinant(a_nan));
  EXPECT_NO_THROW(stan::math::log_determinant(a_infty));
  EXPECT_NO_THROW(stan::math::log_determinant(a_neg_infty));
}
