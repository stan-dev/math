#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/prim.hpp>
#include <test/unit/util.hpp>

TEST(MathMatrixCL, sign_test) {
  Eigen::VectorXd x(7);
  x << 0.2, -0.8, 0, -0, INFINITY, -INFINITY, NAN;

  stan::math::matrix_cl<double> x_cl(x);
  Eigen::VectorXi res = stan::math::from_matrix_cl(stan::math::sign(x_cl));
  EXPECT_MATRIX_EQ(res, stan::math::sign(x));
}

#endif
