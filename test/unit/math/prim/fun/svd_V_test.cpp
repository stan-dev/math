#include <stan/math/prim.hpp>
#include <gtest/gtest.h>
#include <test/unit/util.hpp>

TEST(MathMatrixPrimMat, svd_V) {
    using stan::math::svd_V;
    using stan::math::matrix_d;
    matrix_d m23(2, 3);
    m23 << 1, 3, -5, 7, 9, -11;

    matrix_d m23_V(3, 3);
    m23_V << -0.41176240532160857, 0.81473005032163681, -0.40824829046386291,
      -0.56383954240865775, 0.12417046246885260, 0.81649658092772603,
      0.71591667949570703, 0.56638912538393127, 0.40824829046386307; //Generated using R base::svd
    EXPECT_MATRIX_FLOAT_EQ(m23_V, svd_V(m23));
}
