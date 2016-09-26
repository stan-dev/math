#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>
#include <cstdlib>
#include <ctime>

using Eigen::Matrix;
using Eigen::Dynamic;
using std::rand;

TEST(MathMatrix, matrix_exp_1x1) {
    
    Matrix<double, Dynamic, Dynamic> m1(1,1), m2(1,1);
    m1 << 0;
    m2 << 1;

    expect_matrix_eq(m2, stan::math::matrix_exp(m1));
}

TEST(MathMatrix, matrix_exp_2x2) {

    // example from Moler & Van Loan, 2003
    Matrix<double, Dynamic, Dynamic> m1(2,2), m2(2,2);
    m1 << -49, 24, -64, 31;
    m2 << -.735759, .551819, -1.471518, 1.103638;

    expect_matrix_eq(m2, stan::math::matrix_exp(m1));
}

TEST(MathMatrix, matrix_exp_3x3) {

    // example from http://www.sosmath.com/matrix/expo/expo.html
    Matrix<double, Dynamic, Dynamic> m1(3,3), m2(3,3);
    m1 << 0, 1, 2, 0, 0, -1, 0, 0, 0;
    m2 << 1, 1, 1.5, 0, 1, -1, 0, 0, 1;

    expect_matrix_eq(m2, stan::math::matrix_exp(m1));
}

TEST(MathMatrix, matrix_exp_3x3_2) {

	Matrix<double, Dynamic, Dynamic> m1(3,3), m2(3,3);
	m1 << 89, -66, -18, 20, -14, -4, 360, -270, -73;
	m2 << 245.95891, -182.43047, -49.11821,
		  93.41549, -67.3433, -18.68310,
		  842.54120, -631.90590, -168.14036;

	expect_matrix_eq(m2, stan::math::matrix_exp(m1));
}

TEST(MathMatrix, matrix_exp_100x100) {
	
	int size = 100;
	srand(1); // set seed
	Matrix<double, Dynamic, Dynamic> S = Eigen::MatrixXd::Identity(size, size);
	int col1, col2;
	for(int i = 0; i < 5 * size; i++) {
		col1 = rand() % size;
		col2 = rand() % size;
		while(col1 == col2) col2 = rand() % size;
		S.col(col1) += S.col(col2) * std::pow(-1, rand());
	}
	Matrix<double, 1, Dynamic> diag_elements(size);
	diag_elements.setRandom();
	Matrix<double, 1, Dynamic> exp_diag_elements(size);
	exp_diag_elements = stan::math::exp(diag_elements);
	
	Matrix<double, Dynamic, Dynamic> A;
	A = S * diag_elements.asDiagonal() * stan::math::inverse(S);
	Matrix<double, Dynamic, Dynamic> exp_A, expm_A;
	for(int i=0; i<size; i++) 
	exp_A = S * exp_diag_elements.asDiagonal() * stan::math::inverse(S);
	expm_A = stan::math::matrix_exp(A);
	
	for(int i = 0; i < size; i++)
		for(int j = 0; j < size; j++) {
			if (std::abs(exp_A(i, j)) < 1e-10)
			  EXPECT_NEAR(exp_A(i, j), expm_A(i, j), 1e-2);
			else if (exp_A(i, j) != 0) 
			  EXPECT_FLOAT_EQ(exp_A(i, j), expm_A(i, j));
			else EXPECT_NEAR(exp_A(i, j), expm_A(i, j), 1e-8);
		}
}

TEST(MathMatrix, matrix_exp_exceptions) {

    using stan::math::matrix_exp;

    Matrix<double, Dynamic, Dynamic> m1(0,0), m2(1,2);

    EXPECT_THROW(matrix_exp(m1), std::invalid_argument);
    EXPECT_THROW(matrix_exp(m2), std::invalid_argument);   
}
