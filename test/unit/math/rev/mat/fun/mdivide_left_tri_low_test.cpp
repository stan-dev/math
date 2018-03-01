#include <stan/math/rev/mat.hpp>
#include <stan/math/rev/mat/fun/mdivide_left_tri.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/util.hpp>
#include <test/unit/math/rev/mat/util.hpp>
#include <stdexcept>
#include <iostream>

TEST(AgradRevMatrix, var_var_mdivide_left_tri_low) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::mdivide_left_tri_low;
  using stan::math::vector_v;

  // it only uses a triangular view A

  matrix_v A(2, 2);
  matrix_v B(2, 2);

  A << 1, 0.0 / 0.0, -3, 5;
  B << 2, 5, 12, 109;

  vector_v c(2);
  c << 3, 5;

  matrix_v A_under_B = mdivide_left_tri_low(A, B);
  EXPECT_EQ(2, A_under_B.rows());
  EXPECT_EQ(2, A_under_B.cols());
  for (int i = 0; i < A_under_B.size(); ++i)
    EXPECT_FALSE(stan::math::is_nan(A_under_B(i)));

  EXPECT_NO_THROW(mdivide_left_tri_low(A, B));
  EXPECT_NO_THROW(mdivide_left_tri_low(A, c));

  matrix_v D(3, 2);
  D << 1, 2, 5, -19, 73, 31;
  EXPECT_THROW(mdivide_left_tri_low(D, B), std::invalid_argument);
  EXPECT_THROW(mdivide_left_tri_low(D, c), std::invalid_argument);

  vector_v e(3);
  e << 1, 2, 9;
  EXPECT_THROW(mdivide_left_tri_low(A, e), std::invalid_argument);

  matrix_v A2(2, 2);
  A2 << 1, 0, -3, 5;
  matrix_v B2(2, 2);
  B2 << 2, 5, 12, 109;
  matrix_v A2_under_B2 = A2.inverse() * B2;
  EXPECT_EQ(2, A2_under_B2.rows());
  EXPECT_EQ(2, A2_under_B2.cols());
  for (int i = 0; i < A2_under_B2.size(); ++i)
    EXPECT_FLOAT_EQ(A2_under_B2(i).val(), A_under_B(i).val());
}

TEST(AgradRevMatrix, var_double_mdivide_left_tri_low) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::mdivide_left_tri_low;
  using stan::math::vector_v;

  // it only uses a triangular view A

  matrix_v A(2, 2);
  matrix_d B(2, 2);

  A << 1, 0.0 / 0.0, -3, 5;
  B << 2, 5, 12, 109;

  vector_v c(2);
  c << 3, 5;

  matrix_v A_under_B = mdivide_left_tri_low(A, B);
  EXPECT_EQ(2, A_under_B.rows());
  EXPECT_EQ(2, A_under_B.cols());
  for (int i = 0; i < A_under_B.size(); ++i)
    EXPECT_FALSE(stan::math::is_nan(A_under_B(i)));

  EXPECT_NO_THROW(mdivide_left_tri_low(A, B));
  EXPECT_NO_THROW(mdivide_left_tri_low(A, c));

  matrix_v D(3, 2);
  D << 1, 2, 5, -19, 73, 31;
  EXPECT_THROW(mdivide_left_tri_low(D, B), std::invalid_argument);
  EXPECT_THROW(mdivide_left_tri_low(D, c), std::invalid_argument);

  vector_v e(3);
  e << 1, 2, 9;
  EXPECT_THROW(mdivide_left_tri_low(A, e), std::invalid_argument);

  matrix_v A2(2, 2);
  A2 << 1, 0, -3, 5;
  matrix_v B2(2, 2);
  B2 << 2, 5, 12, 109;
  matrix_v A2_under_B2 = A2.inverse() * B2;
  EXPECT_EQ(2, A2_under_B2.rows());
  EXPECT_EQ(2, A2_under_B2.cols());
  for (int i = 0; i < A2_under_B2.size(); ++i)
    EXPECT_FLOAT_EQ(A2_under_B2(i).val(), A_under_B(i).val());
}

TEST(AgradRevMatrix, double_var_mdivide_left_tri_low) {
  using stan::math::matrix_d;
  using stan::math::matrix_v;
  using stan::math::mdivide_left_tri;
  using stan::math::mdivide_left_tri_low;
  using stan::math::vector_v;

  // it only uses a triangular view A

  matrix_d A(2, 2);
  matrix_v B(2, 2);

  A << 1, 0.0 / 0.0, -3, 5;
  B << 2, 5, 12, 109;

  vector_v c(2);
  c << 3, 5;

  matrix_v A_under_B = mdivide_left_tri_low(A, B);
  EXPECT_EQ(2, A_under_B.rows());
  EXPECT_EQ(2, A_under_B.cols());
  for (int i = 0; i < A_under_B.size(); ++i)
    EXPECT_FALSE(stan::math::is_nan(A_under_B(i)));

  EXPECT_NO_THROW(mdivide_left_tri_low(A, B));
  EXPECT_NO_THROW(mdivide_left_tri_low(A, c));

  matrix_v D(3, 2);
  D << 1, 2, 5, -19, 73, 31;
  EXPECT_THROW(mdivide_left_tri_low(D, B), std::invalid_argument);
  EXPECT_THROW(mdivide_left_tri_low(D, c), std::invalid_argument);

  vector_v e(3);
  e << 1, 2, 9;
  EXPECT_THROW(mdivide_left_tri_low(A, e), std::invalid_argument);

  matrix_v A2(2, 2);
  A2 << 1, 0, -3, 5;
  matrix_v B2(2, 2);
  B2 << 2, 5, 12, 109;
  matrix_v A2_under_B2 = A2.inverse() * B2;
  EXPECT_EQ(2, A2_under_B2.rows());
  EXPECT_EQ(2, A2_under_B2.cols());
  for (int i = 0; i < A2_under_B2.size(); ++i)
    EXPECT_FLOAT_EQ(A2_under_B2(i).val(), A_under_B(i).val());
}
