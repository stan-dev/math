#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>

using Eigen::Matrix;
using Eigen::Dynamic;
using stan::math::expm;

// for dev purposes
using stan::math::factorial;
using stan::math::pow;

TEST(MathMatrix, expm_1x1) {
    
    Matrix<double, Dynamic, Dynamic> m1(1,1), m2(1,1);
    m1 << 0;
    m2 << 1;
        
    expect_matrix_eq(m2, expm(m1));
}

TEST(MathMatrix, expm_2x2) {
    
    // example from Moler & Van Loan, 2003
    Matrix<double, Dynamic, Dynamic> m1(2,2), m2(2,2);
    
    m1 << -49, 24, -64, 31;
    m2 << -.735759, .551819, -1.471518, 1.103638;
    
    expect_matrix_eq(m2, expm(m1));
}

TEST(MathMatrix, expm_symmetric) {

	Matrix<double, Dynamic, Dynamic> m1(3,3), m2(3,3);
	
	m1 << 1, 2, 4, 2, 2, 1, 4, 1, 2;
	
	// result obtained using Pade approximation
	m2 << 263.32996, 177.23225, 273.465,
		  177.23225, 122.44741, 183.11559,
		  273.465, 183.11559, 284.44647;
	
	expect_matrix_eq(m2, expm(m1));

}

TEST(MathMatrix, expm_nilpotent) {

	// example from http://www.sosmath.com/matrix/expo/expo.html
	// For upper-triangular matrix
	Matrix<double, Dynamic, Dynamic> m1(3,3), m2(3,3), m3(3,3), I(3,3);
	m1 << 0, 1, 2, 0, 0, -1, 0, 0, 0;
	m2 << 1, 1, 1.5, 0, 1, -1, 0, 0, 1;  
	m3 << 0, 0, -1, 0, 0, 0, 0, 0, 0;
	I << 1, 0, 0, 0, 1, 0, 0, 0, 1;
	
	EXPECT_EQ(1, factorial(0));
	EXPECT_EQ(1, factorial(1));
	EXPECT_EQ(2, factorial(2));	
	EXPECT_EQ(6, factorial(3));
	
	expect_matrix_eq(I, pow(m1, 0));
	expect_matrix_eq(m1, pow(m1, 1));
	expect_matrix_eq(m3, pow(m1, 2));
	
	// expect_matrix_eq(m2, expm(m1));

}


TEST(MathMatrix, expm_exceptions) {
    
    Matrix<double, Dynamic, Dynamic> m1(0,0), m2(1,2);
    
    EXPECT_THROW(expm(m1), std::invalid_argument);
    EXPECT_THROW(expm(m2), std::invalid_argument);
}