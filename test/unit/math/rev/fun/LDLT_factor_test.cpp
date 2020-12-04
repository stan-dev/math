#include <stan/math/rev/fun/LDLT_factor.hpp>
#include <gtest/gtest.h>

TEST(AgradRevMatrix, LDLT_factor_default_constructor) {
  using stan::math::LDLT_factor;
  using stan::math::var;

  Eigen::Matrix<var, -1, -1> A(2, 2);
  A << 2, 1, 1, 2;

  LDLT_factor<decltype(A)> ldlt_A(A);

  ASSERT_TRUE(ldlt_A.success());
  EXPECT_EQ(2U, ldlt_A.rows());
  EXPECT_EQ(2U, ldlt_A.cols());
  EXPECT_NO_THROW(ldlt_A.solve(A.val()));

  stan::math::recover_memory();
}

TEST(AgradRevMatrix, LDLT_factor_constructor_matrix) {
  using stan::math::LDLT_factor;
  using stan::math::var;

  Eigen::Matrix<var, -1, -1> A(2, 2);
  A << 2, 1, 1, 2;
  Eigen::Matrix<double, -1, -1> B(2, 2);
  B << 2, 1, 1, 2;

  LDLT_factor<decltype(A)> ldlt_A(A);

  ASSERT_TRUE(ldlt_A.success());
  EXPECT_NO_THROW(ldlt_A.vectorD());
  EXPECT_NO_THROW(ldlt_A.solve(B));
  EXPECT_EQ(2U, ldlt_A.rows());
  EXPECT_EQ(2U, ldlt_A.cols());

  stan::math::recover_memory();
}

TEST(AgradRevMatrix, success) {
  using stan::math::LDLT_factor;
  using stan::math::var;

  Eigen::Matrix<var, -1, -1> A(2, 2);
  A << 0, 0, 0, 0;

  LDLT_factor<decltype(A)> ldlt_A(A);
  EXPECT_FALSE(ldlt_A.success());

  A << 2, 1, 1, 2;
  LDLT_factor<decltype(A)> ldlt_A2(A);
  EXPECT_TRUE(ldlt_A2.success());

  stan::math::recover_memory();
}

TEST(AgradRevMatrix, solve) {
  using stan::math::LDLT_factor;
  using stan::math::var;

  Eigen::Matrix<var, -1, -1> A(2, 2);
  A << 2, 1, 1, 2;
  Eigen::Matrix<double, -1, -1> B(2, 2);
  B(0, 0) = 3;
  B(0, 1) = 1;
  B(1, 0) = 2;
  B(1, 1) = 2;

  Eigen::Matrix<double, -1, -1> expected_solve(2, 2);
  expected_solve(0, 0) = 1.333333333;
  expected_solve(0, 1) = 0.0;
  expected_solve(1, 0) = 0.333333333;
  expected_solve(1, 1) = 1.0;
  Eigen::Matrix<double, -1, -1> solve;

  LDLT_factor<decltype(A)> ldlt_A(A);
  ASSERT_TRUE(ldlt_A.success());
  EXPECT_NO_THROW(solve = ldlt_A.solve(B));

  for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++)
      EXPECT_FLOAT_EQ(expected_solve(i, j), solve(i, j));

  stan::math::recover_memory();
}

TEST(AgradRevMatrix, vectorD) {
  using stan::math::LDLT_factor;
  using stan::math::var;

  Eigen::Matrix<var, -1, -1> A(2, 2);
  A << 2, 1, 1, 2;
  Eigen::Matrix<double, -1, -1> vectorD(2, 1);

  Eigen::Matrix<double, -1, -1> A_double(2, 2);
  A_double << 2, 1, 1, 2;
  Eigen::LDLT<Eigen::Matrix<double, 2, 2> > ldlt_double(A_double);
  Eigen::Matrix<double, -1, -1> expected_vectorD(2, 1);
  expected_vectorD = ldlt_double.vectorD();

  LDLT_factor<decltype(A)> ldlt_A(A);
  ASSERT_TRUE(ldlt_A.success());
  EXPECT_NO_THROW(vectorD = ldlt_A.vectorD());
  for (int i = 0; i < 2; i++)
    EXPECT_FLOAT_EQ(expected_vectorD(i), vectorD(i));

  stan::math::recover_memory();
}

TEST(AgradRevMatrix, rows) {
  using stan::math::LDLT_factor;
  using stan::math::var;

  Eigen::Matrix<var, -1, -1> A(2, 2);
  A << 2, 1, 1, 2;

  LDLT_factor<decltype(A)> ldlt_A(A);
  ASSERT_TRUE(ldlt_A.success());
  EXPECT_EQ(2U, ldlt_A.rows());

  stan::math::recover_memory();
}

TEST(AgradRevMatrix, cols) {
  using stan::math::LDLT_factor;
  using stan::math::var;

  Eigen::Matrix<var, -1, -1> A(2, 2);
  A << 2, 1, 1, 2;

  LDLT_factor<decltype(A)> ldlt_A(A);
  ASSERT_TRUE(ldlt_A.success());
  EXPECT_EQ(2U, ldlt_A.cols());

  stan::math::recover_memory();
}

