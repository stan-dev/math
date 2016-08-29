#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>

using stan::math::matrix_fd;
using stan::math::fvar;

// INCOMPLETE - need to finish writing up second test

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
    
    // example from Moler & Van Loan, 2003
    matrix_fd input(2,2), output;
    fvar<double> a, b, c, d;
    
    a.val_ = -49.0;
    b.val_ = 24.0;
    c.val_ = -64.0;
    d.val_ = 31.0;
    a.d_ = 1.0;
    b.d_ = 1.0;
    c.d_ = 1.0;
    d.d_ = 1.0;
    
    input << a, b,
             c, d;
    
    output = matrix_exp(input);
    
    //m2 << -.735759, .551819,
    //-1.471518, 1.103638;
    
    EXPECT_FLOAT_EQ(-0.735759, output(0,0).val_);
    EXPECT_FLOAT_EQ(0.551819, output(0,1).val_);
    EXPECT_FLOAT_EQ(-1.471518, output(1,0).val_);
    EXPECT_FLOAT_EQ(1.103638, output(1,1).val_);
    /*EXPECT_EQ(-0.735759, output(0,0).d_);
    EXPECT_EQ(0.551819, output(0,1).d_);
    EXPECT_EQ(-1.47152, output(1,0).d_);
    EXPECT_EQ(1.103654, output(1,1).d_);*/
    

}
