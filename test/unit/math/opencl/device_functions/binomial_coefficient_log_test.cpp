#ifdef STAN_OPENCL
#include <stan/math.hpp>
#include <stan/math/opencl/kernels/device_functions/lgamma_stirling.hpp>
#include <stan/math/opencl/kernels/device_functions/lgamma_stirling_diff.hpp>
#include <stan/math/opencl/kernels/device_functions/lbeta.hpp>
#include <stan/math/opencl/kernels/device_functions/binomial_coefficient_log.hpp>
#include <stan/math/opencl/kernel_cl.hpp>
#include <gtest/gtest.h>
#include <test/unit/math/expect_near_rel.hpp>
#include <string>

static const std::string test_binomial_coefficient_log_kernel_code
    = STRINGIFY(__kernel void test(__global double *C, __global double *A,
                                   __global double *B) {
        const int i = get_global_id(0);
        C[i] = binomial_coefficient_log(A[i], B[i]);
      });

const stan::math::opencl_kernels::kernel_cl<
    stan::math::opencl_kernels::out_buffer,
    stan::math::opencl_kernels::in_buffer,
    stan::math::opencl_kernels::in_buffer>
    binomial_coefficient_log(
        "test",
        {stan::math::opencl_kernels::lgamma_stirling_device_function,
         stan::math::opencl_kernels::lgamma_stirling_diff_device_function,
         stan::math::opencl_kernels::lbeta_device_function,
         stan::math::opencl_kernels::binomial_coefficient_log_device_function,
         test_binomial_coefficient_log_kernel_code});

TEST(MathMatrixCL, binomial_coefficient_log) {
  Eigen::VectorXd a = Eigen::VectorXd::Random(1000).array() * 5 + 4;
  Eigen::VectorXd b
      = (a.array() + 1) * Eigen::VectorXd::Random(1000).array().abs();

  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> b_cl(b);
  stan::math::matrix_cl<double> res_cl(1000, 1);
  binomial_coefficient_log(cl::NDRange(1000), res_cl, a_cl, b_cl);
  Eigen::VectorXd res = stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);

  EXPECT_NEAR_REL(res, stan::math::binomial_coefficient_log(a, b));
}

TEST(MathMatrixCL, binomial_coefficient_log_edge_cases) {
  Eigen::VectorXd a(5);
  a << NAN, INFINITY, 1.0E50, 5, 5;

  Eigen::VectorXd b(5);
  b << 1, 1, 1, 5, 0;

  stan::math::matrix_cl<double> a_cl(a);
  stan::math::matrix_cl<double> b_cl(b);
  stan::math::matrix_cl<double> res_cl(5, 1);
  binomial_coefficient_log(cl::NDRange(5), res_cl, a_cl, b_cl);
  Eigen::VectorXd res = stan::math::from_matrix_cl<Eigen::VectorXd>(res_cl);

  EXPECT_NEAR_REL(res, stan::math::binomial_coefficient_log(a, b));
}

#endif
