#include <stan/math/fwd/mat.hpp>
#include <gtest/gtest.h>

TEST(MathMatrix, matrix_exp_1x1) {
    
    using stan::math::matrix_fd;
    using stan::math::fvar;
    
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

TEST(MathMatrix, matrix_exp_2x2) {
    
    using stan::math::matrix_fd;
    using stan::math::fvar;
    
    // example from Moler & Van Loan, 2003, section 3
    fvar<double> a, b;
    a.val_ = -1.0;
    b.val_ = -17.0;
    a.d_ = 1.0;
    b.d_ = 1.0;
    
    matrix_fd input(2, 2);
    input << -2*a + 3*b, 1.5*a - 1.5*b,
            -4*a + 4*b, 3*a - 2*b;
    
    matrix_fd output;
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

TEST(MathMatrix, matrix_exp_3x3) {

	using stan::math::matrix_fd;
	using stan::math::fvar;

    fvar<double> a, b, c;
    a.val_ = -1.0;
    b.val_ = 2.0;
    c.val_ = 1.0;
    a.d_ = 1.0;
    b.d_ = 1.0;
    c.d_ = 1.0;
    
    matrix_fd input(3, 3);
    input << -24*a +40*b - 15*c, 18*a - 30*b + 12*c, 5*a - 8*b + 3*c,
    		 20*b - 20*c, -15*b + 16*c, -4*b +4*c,
    		 -120*a + 120*b, 90*a - 90*b, 25*a - 24*b;
    		 
	matrix_fd output;
	output = matrix_exp(input);

	EXPECT_FLOAT_EQ(245.95891, output(0,0).val_);
	EXPECT_FLOAT_EQ(-182.43047, output(0,1).val_);
	EXPECT_FLOAT_EQ(-49.11821, output(0,2).val_);
	EXPECT_FLOAT_EQ(93.41549, output(1,0).val_);
	EXPECT_FLOAT_EQ(-67.3433, output(1,1).val_);
	EXPECT_FLOAT_EQ(-18.68310, output(1,2).val_);
	EXPECT_FLOAT_EQ(842.54120, output(2,0).val_);
	EXPECT_FLOAT_EQ(-631.90590, output(2,1).val_);
	EXPECT_FLOAT_EQ(-168.14036, output(2,2).val_);
	
	// Note: because of way of matrix was constructed, 
	// derivative should be the same as value (same case
	// as in the previous test).  
	
	EXPECT_FLOAT_EQ(245.95891, output(0,0).d_);
	EXPECT_FLOAT_EQ(-182.43047, output(0,1).d_);
	EXPECT_FLOAT_EQ(-49.11821, output(0,2).d_);
	EXPECT_FLOAT_EQ(93.41549, output(1,0).d_);
	EXPECT_FLOAT_EQ(-67.3433, output(1,1).d_);
	EXPECT_FLOAT_EQ(-18.68310, output(1,2).d_);
	EXPECT_FLOAT_EQ(842.54120, output(2,0).d_);
	EXPECT_FLOAT_EQ(-631.90590, output(2,1).d_);
	EXPECT_FLOAT_EQ(-168.14036, output(2,2).d_);
}

TEST(MathMatrix, matrix_exp_exceptions) {
    stan::math::matrix_fd m1(0,0), m2(1,2);
    
    EXPECT_THROW(matrix_exp(m1), std::invalid_argument);
    EXPECT_THROW(matrix_exp(m2), std::invalid_argument);
}
