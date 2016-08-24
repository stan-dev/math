#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/prim/mat/fun/expect_matrix_eq.hpp>

TEST(MathMatrix, matrix_exp) {
    using stan::math::matrix_exp;
    
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> x, y;
    x.resize(1, 1);
    y.resize(1, 1);
    x << 0.0;
    y << 1.0;
        
    expect_matrix_eq(y, matrix_exp(x));
    
    /*
    Eigen::Matrix<stan::math::var, Eigen::Dynamic, Eigen::Dynamic> x_v, y_v;
    x_v.resize(1, 1);
    y_v.resize(1, 1);
    x_v << 0.0;
    y_v << 1.0;
    
    expect_matrix_eq(y_v, matrix_exp(x_v));
    */
    
}
