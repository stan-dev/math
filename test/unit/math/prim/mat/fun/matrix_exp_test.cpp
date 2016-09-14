#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>

TEST(MathMatrix, matrix_exp_1x1) {
    
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(1,1), m2(1,1);
    m1 << 0;
    m2 << 1;
        
    expect_matrix_eq(m2, stan::math::matrix_exp(m1));
}

TEST(MathMatrix, matrix_exp_2x2) {
    
    // example from Moler & Van Loan, 2003
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(2,2), m2(2,2);
    m1 << -49, 24, -64, 31;
    m2 << -.735759, .551819, -1.471518, 1.103638;
    
    expect_matrix_eq(m2, stan::math::matrix_exp(m1));
}

TEST(MathMatrix, matrix_exp_3x3) {
    
    // example from http://www.sosmath.com/matrix/expo/expo.html
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(3,3), m2(3,3);
    m1 << 0, 1, 2, 0, 0, -1, 0, 0, 0;
    m2 << 1, 1, 1.5, 0, 1, -1, 0, 0, 1;
    
    expect_matrix_eq(m2, stan::math::matrix_exp(m1));
}

TEST(MathMatrix, matrix_exp_3x3_2) {

	Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(3,3), m2(3,3);
	m1 << 89, -66, -18, 20, -14, -4, 360, -270, -73;
	m2 << 245.95891, -182.43047, -49.11821,
		  93.41549, -67.3433, -18.68310,
		  842.54120, -631.90590, -168.14036;
		  
	expect_matrix_eq(m2, stan::math::matrix_exp(m1));
}

TEST(MathMatrix, matrix_exp_exceptions) {
    
    using stan::math::matrix_exp;
    
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> m1(0,0), m2(1,2);
    
    EXPECT_THROW(matrix_exp(m1), std::invalid_argument);
    EXPECT_THROW(matrix_exp(m2), std::invalid_argument);
}
