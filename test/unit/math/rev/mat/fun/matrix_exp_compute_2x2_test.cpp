#include <stan/math/rev/mat.hpp>
#include <stan/math/prim/mat/fun/matrix_exp_2x2.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/rev/mat/fun/expect_matrix_eq.hpp>

using stan::math::matrix_v;
using stan::math::matrix_exp_compute_2x2;

TEST(MathMatrix, matrix_exp_compute_2x2_2) {

	matrix_v m1(2,2), m2(2,2);
	
	m1 << 3, 0, 0, 4;
	m2 << exp(3), 0, 0, exp(4);
	
	expect_matrix_eq(m2, matrix_exp_compute_2x2(m1));

}


TEST(MathMatrix, matrix_exp_compute_2x2_1) {
    
    // example from Moler & Van Loan, 2003
    matrix_v m1(2,2), m2(2,2);
    
    m1 << -49, 24, -64, 31;
    m2 << -.735759, .551819, -1.471518, 1.103638;
    
    expect_matrix_eq(m2, matrix_exp_compute_2x2(m1));
}

