#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>

using stan::math::matrix_fd;
using stan::math::fvar;

TEST(MathMatrix, matrix_exp) {
    
    fvar<double> a;
    a.val_ = 0.0;
    a.d_ = 1.0;
    matrix_fd input(1,1), output;
    
    input << a;
    output = matrix_exp(input);
    
    EXPECT_EQ(1.0, output(0,0).val_);
    EXPECT_EQ(1.0, output(0,0).d_);
    
    
    fvar<double> b;
    b.val_ = 1.0;
    b.d_ = 2.0;
    
    input << b;
    output = matrix_exp(input);
    
    EXPECT_EQ(exp(1.0), output(0,0).val_);
    EXPECT_EQ(b.d_ * output(0,0).val_, output(0,0).d_);

}


TEST(MathMatrix, matrix_exp2) {
    
    // example from Moler & Van Loan, 2003, section 3
    matrix_fd input_diag(2,2), input(2,2), output;
    fvar<double> a, b;
    a.val_ = -1.0;
    b.val_ = -17.0;
    a.d_ = 1.0;
    b.d_ = 1.0;
    
    
    input << -2*a + 3*b, 1.5*a - 1.5*b,
            -4*a + 4*b, 3*a - 2*b;
    
    output = matrix_exp(input);
    
    EXPECT_FLOAT_EQ(-0.735759, output(0,0).val_);
    EXPECT_FLOAT_EQ(0.551819, output(0,1).val_);
    EXPECT_FLOAT_EQ(-1.471518, output(1,0).val_);
    EXPECT_FLOAT_EQ(1.103638, output(1,1).val_);
    
    // note: in this particular example, derivatives
    // is the same as the value, due to the way the
    // input matrix was constructed.
    EXPECT_FLOAT_EQ(-0.735759, output(0,0).d_);
    EXPECT_FLOAT_EQ(0.551819, output(0,1).d_);
    EXPECT_FLOAT_EQ(-1.471518, output(1,0).d_);
    EXPECT_FLOAT_EQ(1.103638, output(1,1).d_);

}
