#ifdef STAN_OPENCL
#include <stan/math/prim.hpp>
#include <stan/math/opencl/prim.hpp>
#include <test/unit/util.hpp>

TEST(MathMatrixCL, rank_test) {
  Eigen::VectorXd x(7);
  x << 0.2, -0.8, 0, -0.1, INFINITY, -INFINITY, -1.5;

  stan::math::matrix_cl<double> x_cl(x);
  EXPECT_EQ(stan::math::rank(x_cl, 3), stan::math::rank(x, 3));
  EXPECT_EQ(stan::math::rank(stan::math::transpose(x_cl), 3),
            stan::math::rank(stan::math::transpose(x), 3));
}

#endif
