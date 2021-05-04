#ifdef STAN_OPENCL

#include <stan/math/opencl/kernel_generator.hpp>
#include <stan/math/opencl/matrix_cl.hpp>
#include <stan/math/opencl/copy.hpp>
#include <test/unit/util.hpp>
#include <Eigen/Dense>
#include <gtest/gtest.h>

TEST(KernelGenerator, constant_test) {
  stan::math::matrix_cl<double> m1_cl(stan::math::constant(1.2, 3, 4));

  Eigen::MatrixXd res = stan::math::from_matrix_cl(m1_cl);
  EXPECT_MATRIX_EQ(res, Eigen::MatrixXd::Constant(3, 4, 1.2));
}

#endif
