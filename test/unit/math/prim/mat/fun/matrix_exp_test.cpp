#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>

using Eigen::Matrix;
using Eigen::Dynamic;
using stan::math::matrix_exp;

TEST(MathMatrix, matrix_exp_1) {
    
    Matrix<double, Dynamic, Dynamic> m1(1,1), m2(1,1);
    m1 << 0;
    m2 << 1;
        
    expect_matrix_eq(m2, matrix_exp(m1));
}

TEST(MathMatrix, matrix_exp_2) {
    
    // example from http://www.sosmath.com/matrix/expo/expo.html
    Matrix<double, Dynamic, Dynamic> m1(3,3), m2(3,3), m3(3,3), I(3,3);
    m1 << 0, 1, 2, 0, 0, -1, 0, 0, 0;
    m2 << 1, 1, 1.5, 0, 1, -1, 0, 0, 1;
    
    expect_matrix_eq(m2, matrix_exp(m1));
    
}

TEST(MathMatrix, matrix_exp_2x2) {
    
    // example from Moler & Van Loan, 2003
    Matrix<double, Dynamic, Dynamic> m1(2,2), m2(2,2);
    
    m1 << -49, 24, -64, 31;
    m2 << -.735759, .551819, -1.471518, 1.103638;
    
    expect_matrix_eq(m2, matrix_exp(m1));
}

TEST(MathMatrix, matrix_exp_exceptions) {
    
    Matrix<double, Dynamic, Dynamic> m1(0,0), m2(1,2);
    
    EXPECT_THROW(matrix_exp(m1), std::invalid_argument);
    EXPECT_THROW(matrix_exp(m2), std::invalid_argument);
}