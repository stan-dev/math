#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <stan/math/opencl/kernels/device_functions/log1p_exp.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/expect_near_rel.hpp>
#include <string>

static const std::string test_log1p_exp_kernel_code
    = STRINGIFY(__kernel void test(__global double *B, __global double *A) {
        const int i = get_global_id(0);
        B[i] = log1p_exp(A[i]);
      });

const stan::math::opencl_kernels::kernel_cl<
    stan::math::opencl_kernels::out_buffer,
    stan::math::opencl_kernels::in_buffer>
    log1p_exp("test", {stan::math::opencl_kernels::log1p_exp_device_function,
                       test_log1p_exp_kernel_code});

TEST(MathMatrixCL, log1p_exp) {
  Eigen::VectorXd a = Eigen::VectorXd::Random(1000).array() * 0.9999;

  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> res_cl(1000, 1);
  log1p_exp(cl::NDRange(1000), res_cl, a_cl);
  Eigen::VectorXd res = stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);

  stan::test::expect_near_rel("log1p_exp (OpenCL)", res,
                              stan::math::log1p_exp(a));
}

TEST(MathMatrixCL, log1p_exp_edge_cases) {
  Eigen::VectorXd a(3);
  a << 10000, -10000, NAN;

  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> res_cl(3, 1);
  log1p_exp(cl::NDRange(3), res_cl, a_cl);
  Eigen::VectorXd res = stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);

  stan::test::expect_near_rel("log1p_exp (OpenCL)", res,
                              stan::math::log1p_exp(a));
}

#endif
