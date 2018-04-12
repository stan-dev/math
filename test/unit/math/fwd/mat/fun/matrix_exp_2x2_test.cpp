#include <stan/math/fwd/mat.hpp>
#include <stan/math/prim/mat/fun/matrix_exp_2x2.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, matrix_exp_2x2_1) {
  stan::math::fvar<double> a;
  a.val_ = 3.0;
  a.d_ = 1.0;

  {  // dynamic mat
    stan::math::matrix_fd input(2, 2), output(2, 2);
    input << a, 0, 0, 4;
    output = stan::math::matrix_exp_2x2(input);

    EXPECT_FLOAT_EQ(exp(3), output(0, 0).val_);
    EXPECT_FLOAT_EQ(0, output(1, 0).val_);
    EXPECT_FLOAT_EQ(0, output(0, 1).val_);
    EXPECT_FLOAT_EQ(exp(4), output(1, 1).val_);

    EXPECT_FLOAT_EQ(a.d_ * exp(a.val_), output(0, 0).d_);
    EXPECT_FLOAT_EQ(0, output(1, 0).d_);
    EXPECT_FLOAT_EQ(0, output(0, 1).d_);
    EXPECT_FLOAT_EQ(0, output(1, 1).d_);
  }

  {  // fixed mat
    Eigen::Matrix<stan::math::fvar<double>, 2, 2> input(2, 2), output(2, 2);
    input << a, 0, 0, 4;
    output = stan::math::matrix_exp_2x2(input);

    EXPECT_FLOAT_EQ(exp(3), output(0, 0).val_);
    EXPECT_FLOAT_EQ(0, output(1, 0).val_);
    EXPECT_FLOAT_EQ(0, output(0, 1).val_);
    EXPECT_FLOAT_EQ(exp(4), output(1, 1).val_);

    EXPECT_FLOAT_EQ(a.d_ * exp(a.val_), output(0, 0).d_);
    EXPECT_FLOAT_EQ(0, output(1, 0).d_);
    EXPECT_FLOAT_EQ(0, output(0, 1).d_);
    EXPECT_FLOAT_EQ(0, output(1, 1).d_);
  }
}

TEST(MathMatrix, matrix_exp_2x2_2) {
  // example from Moler & Van Loan, 2003, section 3
  stan::math::fvar<double> a, b;
  a.val_ = -1.0;
  b.val_ = -17.0;
  a.d_ = 1.0;
  b.d_ = 1.0;

  {  // dynamic mat
    stan::math::matrix_fd input(2, 2);
    input << -2 * a + 3 * b, 1.5 * a - 1.5 * b, -4 * a + 4 * b, 3 * a - 2 * b;

    stan::math::matrix_fd output;
    output = stan::math::matrix_exp_2x2(input);

    EXPECT_FLOAT_EQ(-0.735759, output(0, 0).val_);
    EXPECT_FLOAT_EQ(0.551819, output(0, 1).val_);
    EXPECT_FLOAT_EQ(-1.471518, output(1, 0).val_);
    EXPECT_FLOAT_EQ(1.103638, output(1, 1).val_);

    // note: in this particular example, derivatives
    // is the same as the value, due to the way the
    // input matrix was constructed.
    EXPECT_FLOAT_EQ(-0.735759, output(0, 0).d_);
    EXPECT_FLOAT_EQ(0.551819, output(0, 1).d_);
    EXPECT_FLOAT_EQ(-1.471518, output(1, 0).d_);
    EXPECT_FLOAT_EQ(1.103638, output(1, 1).d_);
  }

  {  // fixed mat
    Eigen::Matrix<stan::math::fvar<double>, 2, 2> input;
    input << -2 * a + 3 * b, 1.5 * a - 1.5 * b, -4 * a + 4 * b, 3 * a - 2 * b;

    stan::math::matrix_fd output;
    output = stan::math::matrix_exp_2x2(input);

    EXPECT_FLOAT_EQ(-0.735759, output(0, 0).val_);
    EXPECT_FLOAT_EQ(0.551819, output(0, 1).val_);
    EXPECT_FLOAT_EQ(-1.471518, output(1, 0).val_);
    EXPECT_FLOAT_EQ(1.103638, output(1, 1).val_);

    // note: in this particular example, derivatives
    // is the same as the value, due to the way the
    // input matrix was constructed.
    EXPECT_FLOAT_EQ(-0.735759, output(0, 0).d_);
    EXPECT_FLOAT_EQ(0.551819, output(0, 1).d_);
    EXPECT_FLOAT_EQ(-1.471518, output(1, 0).d_);
    EXPECT_FLOAT_EQ(1.103638, output(1, 1).d_);
  }
}
