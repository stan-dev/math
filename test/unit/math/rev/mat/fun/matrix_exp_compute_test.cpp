#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/mat/fun/MatrixExponential.h>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/expect_matrix_eq.hpp>

using stan::math::matrix_v;
using stan::math::matrix_exp_compute;

TEST(MathMatrix, matrix_exp_compute_1) {
    
    matrix_v m1(1,1), m2(1,1);
    m1 << 0;
    m2 << 1;
        
    expect_matrix_eq(m2, matrix_exp_compute(m1));
}

TEST(MathMatrix, matrix_exp_compute_2) {
    
    // example from Moler & Van Loan, 2003
   matrix_v m1(2,2), m2(2,2);
    
    m1 << -49, 24, -64, 31;
    m2 << -.735759, .551819, -1.471518, 1.103638;
    
    expect_matrix_eq(m2, matrix_exp_compute(m1));
}
