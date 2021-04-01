#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <stan/math/opencl/kernels/device_functions/log1m_exp.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/expect_near_rel.hpp>
#include <string>

static const std::string test_log1m_kernel_code
    = STRINGIFY(__kernel void test(__global double *B, __global double *A) {
        const int i = get_global_id(0);
        B[i] = log1m_exp(A[i]);
      });

const stan::math::opencl_kernels::kernel_cl<
    stan::math::opencl_kernels::out_buffer,
    stan::math::opencl_kernels::in_buffer>
    log1m_exp("test", {stan::math::opencl_kernels::log1m_exp_device_function,
                       test_log1m_kernel_code});

TEST(MathMatrixCL, log1m_exp) {
  Eigen::VectorXd a = Eigen::VectorXd::Random(1000).array() * 0.9999;

  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> res_cl(1000, 1);
  log1m_exp(cl::NDRange(1000), res_cl, a_cl);
  Eigen::VectorXd res = stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);

  stan::test::expect_near_rel("log1m_exp (OpenCL)", res,
                              stan::math::log1m_exp(a));
}

TEST(MathMatrixCL, log1m_exp_edge_cases) {
  Eigen::VectorXd a(11);
  a << -1e10, -1000, -100, -10, -1, -0.1, -1e-5, -1e-10, -1e-20, -1e-40, NAN;

  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> res_cl(11, 1);
  log1m_exp(cl::NDRange(11), res_cl, a_cl);
  Eigen::VectorXd res = stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);

  stan::test::expect_near_rel("log1m_exp (OpenCL)", res,
                              stan::math::log1m_exp(a));
}

#endif
