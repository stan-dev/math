#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <stan/math/opencl/kernels/device_functions/lgamma_stirling.hpp>
#include <stan/math/opencl/kernels/device_functions/lgamma_stirling_diff.hpp>
#include <stan/math/opencl/kernels/device_functions/lbeta.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/expect_near_rel.hpp>
#include <string>

static const std::string test_lbeta_kernel_code = STRINGIFY(__kernel void test(
    __global double *C, __global double *B, __global double *A) {
  const int i = get_global_id(0);
  C[i] = lbeta(A[i], B[i]);
});

const stan::math::opencl_kernels::kernel_cl<
    stan::math::opencl_kernels::out_buffer,
    stan::math::opencl_kernels::in_buffer,
    stan::math::opencl_kernels::in_buffer>
    lbeta("test",
          {stan::math::opencl_kernels::lgamma_stirling_device_function,
           stan::math::opencl_kernels::lgamma_stirling_diff_device_function,
           stan::math::opencl_kernels::lbeta_device_function,
           test_lbeta_kernel_code});

TEST(MathMatrixCL, lbeta) {
  Eigen::VectorXd a = Eigen::VectorXd::Random(1000).array() * 15 + 30;
  Eigen::VectorXd b = Eigen::VectorXd::Random(1000).array() * 15 + 30;

  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> b_cl(b);
  stan::math::matrix_cl<double> res_cl(1000, 1);
  lbeta(cl::NDRange(1000), res_cl, a_cl, b_cl);
  Eigen::VectorXd res = stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);

  EXPECT_NEAR_REL(res, stan::math::lbeta(a, b));
}

TEST(MathMatrixCL, lbeta_edge_cases) {
  Eigen::VectorXd a(3);
  a << NAN, INFINITY, 1.0E50;

  Eigen::VectorXd b(3);
  a << 1, 1, 1;

  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> b_cl(b);
  stan::math::matrix_cl<double> res_cl(3, 1);
  lbeta(cl::NDRange(3), res_cl, a_cl, b_cl);
  Eigen::VectorXd res = stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);

  EXPECT_NEAR_REL(res, stan::math::lbeta(a, b));
}

#endif
