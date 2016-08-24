#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>

// Unit test is incomplete -- still useful to illustrate some bugs

TEST(MathMatrix, matrix_exp) {
    using stan::math::matrix_exp;
    
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x, y;
    x.resize(1, 1);
    y.resize(1, 1);
    x << 0.0;
    y << 1.0;
        
    expect_matrix_eq(y, matrix_exp(x));
    
}
